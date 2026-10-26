

######## Reads table mean, Alpha_diversity & RA calculation ---------------
library(amplicon)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_33sp")
otu_all <- read.csv('field_data_389otu_raw_146.csv', header=T,row.names = 1)
colnames(otu_all)

# OTU data merged by species (average calculated)
otu2 = otu_all[ , -c(1:4, 7:12 )]

otu_by_species <- aggregate(. ~ Plant.species, 
                            data = otu2[, !names(otu2) %in% "Plant.family"], 
                            FUN = mean, na.rm = TRUE)


family_info <- unique(otu2[, c("Plant.species", "Plant.family")])
otu_by_species <- merge(otu_by_species, family_info, by = "Plant.species", all.x = TRUE)


otu_by_species <- otu_by_species[, c("Plant.species", "Plant.family", 
                                     setdiff(names(otu_by_species), 
                                             c("Plant.species", "Plant.family")))]
head(otu_by_species[, 1:5])

data.table::fwrite( otu_by_species , "field_data_389otu_raw_146_33sp.csv",  row.names = TRUE, bom = TRUE)

### 
otu = otu_by_species
rownames(otu) <- otu[, 1]
otu <- otu[ , -c(1,2)]

otu_int <- round(otu)
##2_rarecurve_6_6
rarecurve(otu_int,  step = 10,  cex.lab = 1.2,  main = NULL,
          label = FALSE,  lines = TRUE, col = "tomato", lwd = 0.8,
          xlab = "Number of reads per plant species",  ylab = "Number of OTU" )


#### RA calculation
Tax <- read.csv('field_data_389otu_taxon.csv', header=T, row.names = 1)

####
otu_all <- read.csv('field_data_389otu_raw_146.csv', header=T,row.names = 1)
colnames(otu_all)
rownames(otu_all) <- otu_all[, 12]
otu_t <- t(otu_all[ , -c(1:12)])

format2stamp(otu_t, Tax, thre = 0.0, prefix = "RA_146")



###
otu_all <- read.csv('field_data_389otu_ra_146.csv', header=T,row.names = 1)
colnames(otu_all)

# OTU data merged by species (average calculated)
otu2 = otu_all[ , -c(1:4, 7:12 )]

otu_by_species <- aggregate(. ~ Plant.species, 
                            data = otu2[, !names(otu2) %in% "Plant.family"], 
                            FUN = mean, na.rm = TRUE)


family_info <- unique(otu2[, c("Plant.species", "Plant.family")])
otu_by_species <- merge(otu_by_species, family_info, by = "Plant.species", all.x = TRUE)


otu_by_species <- otu_by_species[, c("Plant.species", "Plant.family", 
                                     setdiff(names(otu_by_species), 
                                             c("Plant.species", "Plant.family")))]
head(otu_by_species[, 1:5])

data.table::fwrite( otu_by_species , "field_data_389otu_ra_146_33sp.csv",  row.names = TRUE, bom = TRUE)


######## Co-occurrence network analysis c -- 146samples raw ---------------------------
#install.packages("glmnet")
#setwd("/Users/zengjing/Downloads") 
#install.packages("SpiecEasi-1.1.1.tar.gz", repos = NULL, type = "source")
# Load required libraries
library(SpiecEasi)
library(phyloseq)
library(igraph)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_33sp")
# Load OTU, taxonomic, and sample data
otu_all <- read.csv("field_data_389otu_raw_146.csv",header=TRUE, row.names=1) 
colnames(otu_all)
rownames(otu_all) <- otu_all[, 12]
otu <- as.matrix( t(otu_all[ , -c(1:12)]) )
tax <- read.csv("field_data_389otu_taxon.csv", header = T, row.names = 1)
sample <- otu_all
rownames(sample) <- sample[, 12]

prev <- rowSums(otu > 0) / ncol(otu)
summary(prev)

table(prev >= 0.05)
table(prev >= 0.10)
table(prev >= 0.15)
table(prev >= 0.20)

c( "all" = sum(prev > 0),
   "5%" = sum(prev >= 0.05),
   "10%" = sum(prev >= 0.10),
   "15%" = sum(prev >= 0.15),
   "20%" = sum(prev >= 0.20)  )

# Create phyloseq object
otu_f <- otu[prev >= 0.05, ]
dim(otu_f)

tax_f <- tax[rownames(otu_f), , drop = FALSE]
all(rownames(otu_f) == rownames(tax_f))

data_f <- phyloseq(
  otu_table(
    otu_f,
    taxa_are_rows = TRUE
  ),
  tax_table(
    as.matrix(tax_f)
  ),
  sample_data(sample)
)

data_f

set.seed(1234)

t1 <- spiec.easi(
  data_f,
  method = "glasso",
  sel.criterion = "stars",
  lambda.log = TRUE,
  lambda.min.ratio = 1e-2,
  nlambda = 20,
  pulsar.select = TRUE,
  pulsar.params = list(
    thresh = 0.05,
    subsample.ratio = 0.8,
    rep.num = 50,
    seed = 1234
  ) )

getOptInd(t1)
opt <- getOptInd(t1)
nlam <- length(t1$lambda)
cat("Optimal lambda index:", opt, "\n")
cat("Number of lambda values:", nlam, "\n")
cat("Optimal lambda:", t1$lambda[opt], "\n")
cat("Minimum lambda:", min(t1$lambda), "\n")

if (opt == nlam) {
  cat("WARNING: Optimal lambda is at the minimum end of the lambda path.\n")
} else {
  cat("Optimal lambda is within the lambda path.\n")
}


adj <- getRefit(t1)

cat(
  "Number of nodes:",
  nrow(adj),
  "\n"
)

cat(
  "Number of edges:",
  Matrix::nnzero(adj) / 2,
  "\n"
)

cat(
  "Network sparsity:",
  sum(adj[upper.tri(adj)] != 0) /
    (nrow(adj) * (nrow(adj) - 1) / 2),
  "\n" )

getOptLambda(t1)
getOptInd(t1)

node_names <- taxa_names(data_f)

g <- graph_from_adjacency_matrix(
  adj, mode = "undirected", diag = FALSE )

# Assign OTU names
V(g)$name <- node_names


g <- delete_vertices(
  g,
  which(igraph::degree(g) == 0)
)

tax_sub <- tax_f[V(g)$name, , drop = FALSE]

# Check that all node names match taxonomy
stopifnot(
  all(V(g)$name %in% rownames(tax_f))
)

V(g)$Kingdom <- tax_sub[V(g)$name, "Kingdom"]
V(g)$Phylum  <- tax_sub[V(g)$name, "Phylum"]
V(g)$Class   <- tax_sub[V(g)$name, "Class"]
V(g)$Order   <- tax_sub[V(g)$name, "Order"]
V(g)$Family  <- tax_sub[V(g)$name, "Family"]
V(g)$Genus   <- tax_sub[V(g)$name, "Genus"]
V(g)$Species <- tax_sub[V(g)$name, "Species"]

adjacency_otu <- as.matrix(
  as_adjacency_matrix(
    g,
    sparse = FALSE
  )
)

write.csv(
  adjacency_otu,
  "1_0916_adjacency_raw_146_spie.csv",
  row.names = TRUE
)

write_graph(
  g,
  file = "1_0916_raw_network_146_spie.graphml",
  format = "graphml"
)

cat(
  "Final network nodes:",
  vcount(g),
  "\n"
)

cat(
  "Final network edges:",
  ecount(g),
  "\n"
)

cat(
  "Network density:",
  edge_density(g),
  "\n"
)


########## 模块提取 
otu_a = read.csv("field_data_389otu_ra_146.csv", header=TRUE, row.names=1) 
rownames(otu_a) <- otu_a[ , 12 ]
otu = otu_a[ , -c(1:12)]
otu_t = t(otu)
head(otu_t)
module0 = read.csv("1_0916_spie_raw_148_447_module.csv", header=TRUE, row.names=1)

# Extract the OTU names from module0
otu_list <- module0$v_name

# Extract the corresponding OTU row from otu_t
otu_module <- otu_t[rownames(otu_t) %in% otu_list, ]

write.csv( otu_module , "1_0916_spie_raw_148_447_module0.csv", row.names = T)

####
module = read.csv("1_0916_spie_raw_148_447_module1.csv", header=TRUE, row.names=1)

module_by_species <- aggregate(. ~ modularity_class, data = module , FUN = sum , na.rm = TRUE)

head(module_by_species[, 1:5])

write.csv(module_by_species, "1_0916_spie_raw_148_447_module2.csv", row.names = FALSE)



######## --mean cal-- 146 33species ----------------
library(plyr) 
library(ggtree)
library(ape) 
library(picante) 
library(phylosignal) 
library(phylobase) 
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_33sp") 
data = read.csv("field_data_146.csv",header=TRUE, row.names=1)

data_mean <- ddply(data, c("Plant.species","Plant.species0",  "Plant.family"),summarise,   
                   P.concentrations = mean(P.concentrations,na.rm=TRUE),
                   Root.diameter = mean(Root.diameter,na.rm=TRUE),
                   Root.length.colonization = mean(Root.length.colonization,na.rm=TRUE),
                   
                   Glomeraceae = mean(Glomeraceae,na.rm=TRUE),
                   Acaulosporaceae = mean(Acaulosporaceae,na.rm=TRUE),
                   
                   Dominikia = mean(Dominikia, na.rm=TRUE),
                   Rhizophagus = mean(Rhizophagus,na.rm=TRUE),
                   Acaulospora = mean( Acaulospora, na.rm=TRUE),
                   Polonospora = mean( Polonospora , na.rm=TRUE),
                   Funneliformis = mean( Funneliformis, na.rm=TRUE),
                   Entrophospora = mean( Entrophospora , na.rm=TRUE),
                   Kamienskia = mean( Kamienskia , na.rm=TRUE),
                   Diversispora = mean( Diversispora , na.rm=TRUE),
                   Others = mean( Others , na.rm=TRUE),
                   
                   OTU1 = mean(OTU1, na.rm=TRUE),
                   OTU2 = mean(OTU2, na.rm=TRUE),
                   
                   
                   Family.richness = mean(Family.richness,na.rm=TRUE),
                   Genus.richness = mean(Genus.richness,na.rm=TRUE),
                   OTU.richness = mean(OTU.richness, na.rm=TRUE),
                   Module.1 = mean(Module.1, na.rm=TRUE),
                   Module.2 = mean(Module.2, na.rm=TRUE),
                   Module.3 = mean( Module.3 , na.rm=TRUE),
                   Module.4 = mean( Module.4 , na.rm=TRUE),
                   Module.5 = mean( Module.5 , na.rm=TRUE)
                   
)
data_mean


data.table::fwrite(data_mean,"field_data_146_33sp_mean.csv",  row.names = TRUE, bom = TRUE)


########  v.phylomaker tree  33sp  --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_35sp")
df <- read.csv("33sp.csv")
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

write.tree(sp_tr, "1_33sp_phylomakers1.nwk")

#sp_tr$node.label <- NULL
#write.tree(sp_tr, "1-35sp-phylomakers1-mega.nwk")

is.rooted(sp_tr)       
is.binary.tree(sp_tr)  

####
library(dplyr)
tree <- read.tree("1_33sp_phylomakers1.nwk")
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

write.tree( tree, file = "1_33sp_phylomakers1_abbr.nwk")


#########################################################################################################
######## -- Fig 2 -- Relationship between RD and AM fungal variables ----------------------
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
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_33sp")
data = read.csv("field_data_146_33sp.csv",header=TRUE, row.names=1) 
colnames(data)

##
rooted_tree <- read.tree("1_33sp_phylomakers1_abbr.nwk") 
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

fig2a <- ggplot( data , aes(y = Rhizophagus , x = Root.diameter  ))+   labs(title = "Field: 33 species") +     # plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0)
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = expression("RA of "*italic(Rhizophagus)*" (%)")) +
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22",  show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth= 0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci1,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model11 ), linewidth = 0.8 , color = "#CB4641", linetype = "solid" ) + 
  theme_classic() + theme(plot.title = element_blank(),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "grey22", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 1.01, color = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.269, 0.002 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 2.1, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.269, 0.002 ), size = 4, inherit.aes = FALSE )
fig2a



##
lm_model <- lm( OTU.richness ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  OTU.richness ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model2 <- predict(pgls_model, newdata = data )  

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
pgls_ci2 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci2)

fig2b <- ggplot( data , aes(y = OTU.richness , x = Root.diameter  ))+      labs(title = "Field: 33 species") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "OTU richness" )+ 
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid"  )+ 
  #geom_ribbon(data = pgls_ci2 ,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model2 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_blank() ,axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 1.01, colour = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.051, 0.207 ), size = 4 , inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 2.1, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.051, 0.207 ), size = 4 , inherit.aes = FALSE )
fig2b


##
lm_model <- lm( Family.richness ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  Family.richness ~ Root.diameter , data = comp_data, lambda = "ML") 
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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( Family.richness ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci3 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci3)

fig2c <- ggplot( data , aes(y = Family.richness , x = Root.diameter  ))+      labs(title = "Field: 33 species") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "Family richness" ) +
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid"  )+ 
  #geom_ribbon(data = pgls_ci3,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model3 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_blank(),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = 1.01 , colour = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.016, 0.490 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = 2.1 , color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.016, 0.490 ), size = 4, inherit.aes = FALSE )
fig2c


##
lm_model <- lm( PCoA1 ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  PCoA1 ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model4 <- predict(pgls_model, newdata = data )  


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
pgls_ci4 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci4)

fig2d <- ggplot( data , aes( y = PCoA1 , x = Root.diameter  ))+      labs(title = "Field: 33 species") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "PCoA axis1 of community composition" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22",  show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "black",  fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "dashed"  )+ 
  #geom_ribbon(data = pgls_ci4,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model4 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_blank() ,axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust =  3.0 , colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.048, 0.222 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust =  4.0 , color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.011, 0.561 ), size = 4, inherit.aes = FALSE )
fig2d


##
lm_model <- lm( PCoA2 ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  PCoA2 ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model5 <- predict(pgls_model, newdata = data )  


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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( PCoA2 ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci5 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci5)

fig2e <- ggplot( data , aes( y = PCoA2 , x = Root.diameter  ))+      labs(title = "Field: 33 species") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "PCoA axis2 of community composition" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22",   show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black",  fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci5,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model5 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_blank() ,axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust =  1.0 , colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.309, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust =  2.0, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.309, 0.001 ), size = 4, inherit.aes = FALSE )
fig2e


##
lm_model <- lm( Module.2 ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  Module.2 ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model6 <- predict(pgls_model, newdata = data )  

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
pgls_ci6 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci6)

fig2f <- ggplot( data , aes( y = Module.2 , x = Root.diameter ))+       labs(title = "Field: 33 species") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "RA of Module #2 (%)" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22",   show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci6,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model6 ), linewidth = 0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_blank(), axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 2.2, colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.267, 0.002 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 3.2, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.267, 0.002 ), size = 4, inherit.aes = FALSE )
fig2f


#
fig2 <- cowplot::plot_grid( fig2a,  fig2b, fig2c, 
                            fig2d, fig2e, fig2f,
                            
                            ncol = 3, nrow = 2, labels = c( "(a)", "(b)",'(c)',
                                                            "(d)", "(e)","(f)",  
                                                            "(g)", "(h)", "(i)",
                                                            "(j)", "(k)", "(l)" ), 
                            label_fontface = "bold", label_x = -0.02, label_y = 1.04 , align = c("hv"))
fig2

setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom")
eoffice::topptx( fig2 , width = 9 ,  height = 7 , "1_fig2_146_33sp_9_7.pptx" ) 


######## -- Fig 3 -- Relationship between plant traits and AM fungal variables ----------------------
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
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_33sp")
data = read.csv("field_data_146_33sp.csv",header=TRUE, row.names=1)  
colnames(data)

##
rooted_tree <- read.tree("1_33sp_phylomakers1_abbr.nwk") 
plot(rooted_tree)

rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data, names.col = "Plant.species", vcv = T, vcv.dim = 2 )
nrow(comp_data$data)
nrow(data)
setdiff(data$Plant.species, rownames(comp_data$data))


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


fig3a <- ggplot( data, aes(y = P.concentrations , x = Rhizophagus  ))+    labs(title = "Field: 33 species") + # plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0)
  scale_x_continuous(name = expression("RA of "*italic(Rhizophagus)*" (%)")) + scale_y_continuous(name = expression("Shoot P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22",  show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid" )+ 
  #geom_ribbon(data = pgls_ci2,  aes(x = Rhizophagus, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model22 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_blank() ,axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  #guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 1.01, colour = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.131, 0.038), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 2.1, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.070, 0.137), size = 4, inherit.aes = FALSE )

fig3a


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

fig3b <- ggplot( data, aes(y = P.concentrations , x = OTU.richness  ))+     labs(title = "Field: 33 species") +  
  scale_x_continuous(name ="OTU richness") + scale_y_continuous(name = expression("Shoot P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,   shape = 21,  fill = "#d3deea" , color = "grey22",  show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid" )+ 
  #geom_ribbon(data = pgls_ci4,  aes(x = OTU.richness, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model4 ), linewidth = 0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_blank() ,axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 1.01, colour = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.001, 0.877 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 2.1, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.002, 0.797 ), size = 4, inherit.aes = FALSE )

fig3b



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


fig3c <- ggplot( data, aes(y = P.concentrations ,  x = Family.richness  )) +    labs(title = "Field: 33 species") +   
  scale_x_continuous(name ="Family richness") + scale_y_continuous(name = expression("Shoot P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,   shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci6,  aes(x = Family.richness, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model6 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_blank(), axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust =  1.01, colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.138, 0.033 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust =  2.1, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.132, 0.037 ), size = 4, inherit.aes = FALSE )

fig3c


###
lm_model <- lm( P.concentrations ~ PCoA1, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ PCoA1 , data = comp_data, lambda = "ML") 
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
pred_grid <- seq(min(dat0$PCoA1), max(dat0$PCoA1), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( P.concentrations ~ PCoA1 ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( PCoA1 = pred_grid) ) } } }

# CI
pgls_ci8 <- data.frame( PCoA1 = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig3d <- ggplot( data, aes( y = P.concentrations , x = PCoA1  )) +     labs(title = "Field: 33 species") +
  scale_x_continuous(name ="PCoA axis1 of community composition") + scale_y_continuous(name = expression("Shoot P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,   shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "black",  fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "dashed" )+ 
  geom_ribbon(data = pgls_ci8,  aes(x = PCoA1, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model8 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_blank(),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text = element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust =  1.01, colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.070, 0.136 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust =  2.1, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.202, 0.009 ), size = 4, inherit.aes = FALSE )


fig3d


###
lm_model <- lm( P.concentrations ~ PCoA2, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ PCoA2 , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model9 <- predict(pgls_model, newdata = data )  

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
pgls_ci9 <- data.frame( PCoA2 = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig3e <- ggplot( data, aes(y = P.concentrations , x = PCoA2  )) +     labs(title = "Field: 33 species") +  
  scale_x_continuous(name ="PCoA axis2 of community composition") + scale_y_continuous(name = expression("Shoot P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,   shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black",  fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci9,  aes(x = PCoA2, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model9 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_blank(),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text = element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust =  1.01, colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.194, 0.010 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust =  2.1, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.130, 0.039 ), size = 4, inherit.aes = FALSE )


fig3e


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


fig3f <- ggplot( data, aes(y = P.concentrations , x = Module.2  )) +   labs(title = "Field: 33 species") +    
  scale_x_continuous(name ="RA of Module #2 (%)") + scale_y_continuous(name = expression("Shoot P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,   shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci10,  aes(x = Module.2, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model222 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_blank(), axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = 1.01, colour = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.333, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = 2.1, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.280, 0.002 ), size = 4, inherit.aes = FALSE )

fig3f


fig3 <- cowplot::plot_grid( fig3a, fig3f, fig3b, fig3c, 
                            fig3d, fig3e, 
                            
                            ncol = 3, nrow = 2, labels = c( "(a)", "(b)",'(c)',
                                                            "(d)", "(e)","(f)",  
                                                            "(g)", "(h)", "(i)",
                                                            "(j)", "(k)", "(l)" ), 
                            label_fontface = "bold", label_x = -0.02, label_y = 1.01 , align = c("hv"))
fig3

setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom")
eoffice::topptx( fig3 , width = 9 ,  height = 7 , "1_fig3_146_33sp_9_7.pptx" ) 


#########################################################################################################
######## -- Fig S3 a -- AM fungal genus ---------------
library(ggalluvial)
library(ggplot2)
library(tidyr)
library(dplyr)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_33sp")
data = read.csv("field_data_146_33sp.csv",header=TRUE, row.names=1) 
colnames(data)

genus = data[, c(9:17)]
head(genus)
rownames(genus) <- data[, 1]
genus$Sample <- rownames(genus)

genus_long <- genus  %>%  pivot_longer(cols = -Sample, names_to = "Genus", values_to = "Abundance")

genus_order <- c("Dominikia", "Rhizophagus", "Acaulospora", "Polonospora", 
                 "Funneliformis", "Entrophospora",   "Kamienskia", "Diversispora", "Others")

genus_long$Genus <- factor(genus_long$Genus, levels = genus_order)

genus_labels <- setNames(
  c(expression(italic("Dominikia")), 
    expression(italic("Rhizophagus")), 
    expression(italic("Acaulospora")), 
    expression(italic("Polonospora")), 
    expression(italic("Funneliformis")), 
    expression(italic("Entrophospora")), 
    expression(italic("Kamienskia")), 
    expression(italic("Diversispora")), 
    "Others"),    c("Dominikia", "Rhizophagus", "Acaulospora", "Polonospora", 
                    "Funneliformis", "Entrophospora",   "Kamienskia", "Diversispora", "Others") )


p_sankey1 <- ggplot(genus_long, aes(x = Sample, y = Abundance, alluvium = Genus,  stratum = Genus, fill = Genus) ) +
  geom_alluvium(aes(fill = Genus),   alpha = 0.3,  width = 0.5, show.legend = F ) +
  geom_stratum(width = 0.5, alpha = 0.7, color = "black",show.legend = T ) +
  scale_fill_manual(values = c("Dominikia" = "#5880ba",  "Rhizophagus" = "#E9746A", "Acaulospora" = "#7B8F5E",
                               "Polonospora" = "#B07AA8","Funneliformis" = "#497d75", "Entrophospora" = "#FFA64D",
                               "Kamienskia" = "#A65928",   "Diversispora" = "#8FB9D1","Others" = "#999999"), 
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


######## -- Fig S3 b -- pcoa ------------------------------------
library(vegan)
library(ggplot2)
library(ggrepel)
library(ape)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_33sp") 
data = read.csv("field_data_146_33sp.csv",header=TRUE, row.names=1) 
otu_a <- read.csv("field_data_389otu_raw_146_33sp.csv",header=TRUE, row.names=1)
otu = otu_a[, -c(1,2)] 
rownames(otu) = otu_a[, 1]

vare.dist <- vegdist(otu, method="bray")

pcoa_res <- pcoa(vare.dist)

pcoa_axis <- as.data.frame(pcoa_res$vectors[, 1:2])

colnames(pcoa_axis) <- c("PCoA1", "PCoA2")

eig_percent <- round(pcoa_res$values$Relative_eig[1:2] * 100, 2)
eig_percent

write.csv( pcoa_axis , '1_0916_raw_146_33sp_pcoa12.csv')  

###
set.seed(123)
# envfit
fit <- envfit( pcoa_res$vectors , data$Root.diameter, permutations = 9999)
fit
fit <- envfit( pcoa_res$vectors , data$P.concentrations, permutations = 9999)
fit

adonis2(vare.dist ~ data$Root.diameter, permutations = 9999) 
adonis2(vare.dist ~ data$P.concentrations, permutations = 9999) 
#  Mantel test
rd_dist <- dist(data$P.concentrations, method = "euclidean")
set.seed(123)
mantel_rd <- mantel(  vare.dist,  rd_dist,  method = "pearson",  permutations = 9999 )
mantel_rd

fit <- envfit( pcoa_res$vectors, data[, c("Root.diameter", "P.concentrations")], permutations = 9999 ) 

# 
fit_df <- data.frame( Axis1 = fit$vectors$arrows[, 1], Axis2 = fit$vectors$arrows[, 2],
                      r2 = fit$vectors$r,
                      P = fit$vectors$pvals )

fit_df$Variable <- rownames(fit$vectors$arrows)

# 调整箭头长度（可根据图形再微调）
arrow_mul <- 0.6
fit_df$PCoA1 <- fit_df$Axis1 * arrow_mul
fit_df$PCoA2 <- fit_df$Axis2 * arrow_mul
fit_df$Variable <- c("Root diameter", "Tissue P concentration")

# 准备数据
pcoa_axis$Sample <- rownames(otu)
eig_percent_formatted <- sprintf("%.2f", eig_percent)

# 基础 PCoA 图
pcoa <- ggplot(pcoa_axis, aes(x = PCoA1, y = PCoA2, )) + 
  geom_point( size = 2.5, shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
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


######## -- Fig S3 c -- module genus ----------------------------------
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(tidyr)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_33sp")
data1 = read.csv("field_data_module5.csv",header=TRUE, row.names=1) 
head(data1)

Module = c("Module #1", "Module #2", "Module #3", "Module #4", "Module #5")
data1$Module <- Module 
Module <- rownames(data1) 

data_long <- data1 %>%
  pivot_longer(cols = -Module, names_to = "Genus", values_to = "Abundance")

# 
genus_order <- c("Dominikia", "Rhizophagus", "Acaulospora", "Polonospora", 
                 "Funneliformis", "Entrophospora", "Kamienskia",  "Diversispora", "Others")

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
    expression(italic("Diversispora")), 
    "Others"),  # Others 不加 italic
  c("Dominikia", "Rhizophagus", "Acaulospora", "Polonospora", 
    "Funneliformis", "Entrophospora", "Kamienskia", "Diversispora", "Others") )


p_sankey2 <- ggplot(data_long, aes(x = Module, y = Abundance, alluvium = Genus,  stratum = Genus, fill = Genus) ) +
  geom_alluvium(aes(fill = Genus),   alpha = 0.3,  width = 0.3, show.legend = F ) +
  geom_stratum(width = 0.3, alpha = 0.7, show.legend = T ) +
  scale_fill_manual(values = c("Dominikia" = "#5880ba",  "Rhizophagus" = "#E9746A", "Acaulospora" = "#7B8F5E",
                               "Polonospora" = "#B07AA8","Funneliformis" = "#497d75", "Entrophospora" = "#FFA64D",
                               "Kamienskia" = "#A65928",   "Diversispora" = "#8FB9D1","Others" = "#999999"), 
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

setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom")
eoffice::topptx( figs5.2 , width = 14 ,  height= 8.5 , "1_figs3_146_33sp_14_8.5.pptx" )


######## -- Fig S4 -- Distribution and sharing patterns of AM fungal OTUs ---------------
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_33sp")
data_all = read.csv("field_data_389otu_raw_146_33sp.csv", header=TRUE, row.names=1) 
colnames(data_all)
rownames(data_all) <- data_all[, 1]
data_all <- data_all[ , -c(1,2)]
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
tax = read.csv("field_data_389otu_taxon.csv", header= T, row.names=1 ) 
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
length(unique(sankey_data$Plant)) # 27species

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
  labs( x = NULL, y = "Interaction frequency", fill = "Genus", title = "Distribution of 88 plant-specific AM fungal OTUs" ) +
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

setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom")
eoffice::topptx( figs4 , width = 12 ,  height = 16 , "1_figs4_146_33sp_8_20.pptx" ) 


######## -- Table S2 -- Plant phylogenetic signals of AM fungal variables 33species ----------------
library(plyr) 
library(ggtree)
library(ape) 
library(picante) 
library(phylosignal) 
library(phylobase) 
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_33sp")  
data_mean = read.csv("field_data_146_33sp.csv",header=TRUE, row.names=1) 

data_mean1 = data_mean
rownames(data_mean1) <- data_mean[, 1]
data_mean1 <- data_mean1[ , -c(1,2,3)]

##
rooted_tree <- read.tree("1_33sp_phylomakers1_abbr.nwk") 
plot(rooted_tree)

tree_p <- as(rooted_tree, "phylo4")
class(tree_p)
p4d <- phylo4d(tree_p, data_mean1)
phy.sig <- phyloSignal(p4d, reps = 9999)
phy.sig ## 35sp 150

write.csv(phy.sig, '1_1024_146_33sp_physig.csv')  



######## Relationship between plant phylogeny and AM fungal community composition ---------------
library(vegan)
library(picante)
library(ggtree)
library(ape)
library(phytools)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_33sp")
data = read.csv("field_data_146_33sp.csv",header=TRUE, row.names=1) 
otu_a <- read.csv("field_data_389otu_raw_146_33sp.csv",header=TRUE, row.names=1)
otu = otu_a[, -c(1,2)] 
rownames(otu) = otu_a[, 1]

set.seed(123)
vare.dist <- vegdist(otu, method="bray")

#adonis2(vare.dist ~ data$Plant.species, permutations = 999) 

# Relationship between plant phylogeny and AM fungal communities
phylo <- read.tree("1_33sp_phylomakers1_abbr.nwk")
plot(phylo)

phylo_dist <- cophenetic.phylo(phylo)  

sample_species <- data$Plant.species
sample_phylo_dist <- phylo_dist[sample_species, sample_species]
all(sample_species %in% rownames(phylo_dist))  
dim(sample_phylo_dist)  

head(rownames(sample_phylo_dist))
head(sample_species)
all(rownames(sample_phylo_dist) == sample_species)  

matrix_data <- as.matrix(sample_phylo_dist)
dist.phy <- as.dist(matrix_data)
otu_phy <- mantel(vare.dist, dist.phy,  method = 'pearson', permutations = 999, na.rm = T)
otu_phy

# Relationship between plant traits and AM fungal communities
set.seed(123)
vare.mds <- monoMDS(vare.dist, model = c("local"), threshold = 0.8,  k=2, 
                    autotransform=T, stress = 1, maxit = 200, sratmax = 0.9999, parallel=2,  noshare=F) 

vare.mds

plot(vare.mds, type = "n") 
points(vare.mds,  pch = 19, col = "grey")
fit1 <- ordisurf(vare.mds ~ Root.diameter, data , col = "blue", add = T, select = F, method = "GACV.Cp")  #  REML
summary(fit1) 

plot(vare.mds, type = "n") 
points(vare.mds,  pch = 19, col = "grey")
fit2 <- ordisurf(vare.mds ~ P.concentrations, data , col = "blue", add = T, select = F, method = "GACV.Cp")  #  REML
summary(fit2) 


######## -- Nestedness and modularity of plant-AM fungal bipartite network ---------------
library(bipartite) 
library(vegan)
library(dplyr) 
rm(list = ls()) 
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_33sp")  
# nestedness of plant-AM fungal bipartite network
otu_a = read.csv("field_data_389otu_raw_146_33sp.csv",header=TRUE, row.names=1) 

binary_matrix <- ifelse(otu_a > 0, 1, 0)

set.seed(1234)
null_model <- nullmodel( binary_matrix , N = 1000, method = "r2dtable")

nested_NODF <- function(x) nested(x, method = "NODF")

result <- oecosimu( binary_matrix, nested_NODF, method = "r2dtable", nsimul = 1000 )
result

result$oecosimu$statistic
#statistic    SES   mean   2.5%    50%  97.5% Pr(sim.)    
#NODF    31.961 23.748 23.482 22.792 23.482 24.188 0.000999 ***
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

result$oecosimu$means
#  23.48156
sd(result$oecosimu$simulated)
#  0.357062


# modularity of plant-AM fungal bipartite network
set.seed(1234)
otu_a = read.csv("field_data_389otu_raw_146_33sp.csv",header=TRUE, row.names=1) 

binary_matrix <- ifelse(otu_a > 0, 1, 0)

mod_result <- computeModules(binary_matrix)
plotModuleWeb(mod_result)

mod_obs <- mod_result
mod_obs@likelihood
Q_obs <- mod_obs@likelihood
Q_obs
# [1] 0.2077317

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
# [1] 0.2214986
sd(Q_null)
# [1] 0.003022502

SES <- (Q_obs - mean(Q_null))/sd(Q_null)
SES
# [1] -4.554827
p_two <- (sum(abs(Q_null - mean(Q_null)) >= abs(Q_obs - mean(Q_null))) + 1)/(n_sim + 1)
p_two
# [1] 0.000999001
p_upper <- (sum(Q_null >= Q_obs) + 1 )/(n_sim+1)
p_upper
# [1] 1



######## -- R2 P ols  pgls  33sp ---------------
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
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_33sp")
data = read.csv("field_data_146_33sp.csv",header=TRUE, row.names=1) 
colnames(data)

##
rooted_tree <- read.tree("1_33sp_phylomakers1_abbr.nwk") 
plot(rooted_tree)

rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data, names.col = "Plant.species", vcv = T, vcv.dim = 2 )
nrow(comp_data$data)
nrow(data)
setdiff(data$Plant.species, rownames(comp_data$data))


##
lm_model <- lm( Polonospora ~ Root.diameter, data = data ) 
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


#########################################################################################################
######### Sensitivity analysis1 28 sp ###################################################################
######## Statistics on community sequence -------------------------
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_28sp")
otu_all = read.csv("field_data_391otu_raw_155.csv", header=TRUE, row.names=1)
colnames(otu_all)
otu = otu_all[ , -c( 1:14 )]

read_counts <- rowSums(otu)
summary(read_counts)
hist(read_counts, breaks = 100, main = "Sequencing depth distribution")

five_quantiles <- quantile(read_counts, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1))
print(five_quantiles)

# Mark quantiles on the histogram.
hist(read_counts, breaks = 20 , main = "Sequencing depth distribution")
abline(v = five_quantiles, col = "red", lty = 2, lwd =  0.8 )
legend("topright", legend = names(five_quantiles), col = "red", lty = 2, lwd = 1.5, cex = 0.7)

# Calculate the total number of sequences for each sample.
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


######## -- rrarefy & rarecurve plot —— 118samples 28 sp  ---------------
library(vegan)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_28sp")
otu_all = read.csv("field_data_386otu_raw_118.csv", header=TRUE, row.names=1)
colnames(otu_all)
otu1 = otu_all[ , -c( 1:12 )]

par(mar = c(5, 5, 1, 2))
rarecurve( otu1 ,  step = 10,  cex.lab = 1.2,  main = "122samples" ,
           label = F ,  lines = TRUE, col = "steelblue", lwd = 0.8,
           xlab = "Number of reads per sample",  ylab = "Number of OTU" )

###
set.seed(123)
otu_rrarefy = rrarefy(otu1, min(rowSums(otu1)))   # 9908

otu_all_rrarefy = cbind(otu_all[1:12], otu_rrarefy)

#data.table::fwrite(otu_all_rrarefy, "field_data_386otu_rare_118.csv",  row.names = TRUE, bom = TRUE)

##2_rarecurve_122samples_6_6
par(mar = c(5, 5, 1, 2))
rarecurve(otu_rrarefy,  step = 10,  cex.lab = 1.2,  main = NULL,
          label = FALSE,  lines = TRUE, col = "steelblue", lwd = 0.8,
          xlab = "Number of reads per sample",  ylab = "Number of OTU" )


# 
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

#data.table::fwrite( otu_by_species , "field_data_386otu_rare_118_28sp.csv",  row.names = TRUE, bom = TRUE)


otu = otu_by_species
rownames(otu) <- otu[, 1]
otu <- otu[ , -c(1,2)]

otu_int <- round(otu)
##2_rarecurve_118_28species_6_6
rarecurve(otu_int,  step = 10,  cex.lab = 1.2,  main = NULL,
          label = FALSE,  lines = TRUE, col = "tomato", lwd = 0.8,
          xlab = "Number of reads per plant species",  ylab = "Number of OTU" )



######## CA Reads table mean  & Alpha_diversity & RA calculation ---------------
library(amplicon)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_28sp")
otu_all <- read.csv('field_data_386otu_rare_118.csv', header=T,row.names = 1)
colnames(otu_all)

# OTU data merged by species (average calculated)
otu2 = otu_all[ , -c(1:4, 7:12 )]

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

### 
otu = otu_by_species
rownames(otu) <- otu[, 1]
otu <- otu[ , -c(1,2)]

otu_int <- round(otu)
##2_rarecurve_6_6
rarecurve(otu_int,  step = 10,  cex.lab = 1.2,  main = NULL,
          label = FALSE,  lines = TRUE, col = "tomato", lwd = 0.8,
          xlab = "Number of reads per plant species",  ylab = "Number of OTU" )


#### RA calculation
##
Tax <- read.csv('field_data_389otu_taxon.csv', header=T, row.names = 1)

####
####
otu_all <- read.csv('field_data_386otu_rare_118.csv', header=T,row.names = 1)
colnames(otu_all)
rownames(otu_all) <- otu_all[, 12]
otu_t <- t(otu_all[ , -c(1:12)])

format2stamp(otu_t, Tax, thre = 0.0, prefix = "RA_118")



######## CA Co-occurrence network analysis SPIEC-EASI -- 118samples rare ---------------------------
#install.packages("glmnet")
#setwd("/Users/zengjing/Downloads") 
#install.packages("SpiecEasi-1.1.1.tar.gz", repos = NULL, type = "source")
# Load required libraries
library(SpiecEasi)
library(phyloseq)
library(igraph)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_28sp")
# Load OTU, taxonomic, and sample data
otu_all <- read.csv("field_data_386otu_rare_118.csv",header=TRUE, row.names=1) 
rownames(otu_all) <- otu_all[, 12]
otu <- as.matrix( t(otu_all[,-c(1:12)]) )
tax <- read.csv("field_data_389otu_taxon.csv", header = T, row.names = 1)
sample <- otu_all
rownames(sample) <- sample[, 12]

prev <- rowSums(otu > 0) / ncol(otu)
summary(prev)

table(prev >= 0.05)
table(prev >= 0.10)
table(prev >= 0.15)
table(prev >= 0.20)

c( "all" = sum(prev > 0),
   "5%" = sum(prev >= 0.05),
   "10%" = sum(prev >= 0.10),
   "15%" = sum(prev >= 0.15),
   "20%" = sum(prev >= 0.20)  )

# Create phyloseq object
otu_f <- otu[prev >= 0.05, ]
dim(otu_f)

tax_f <- tax[rownames(otu_f), , drop = FALSE]
all(rownames(otu_f) == rownames(tax_f))

data_f <- phyloseq(
  otu_table(
    otu_f,
    taxa_are_rows = TRUE
  ),
  tax_table(
    as.matrix(tax_f)
  ),
  sample_data(sample)
)

data_f

set.seed(1234)
t1 <- spiec.easi(
  data_f,
  method = "glasso",
  sel.criterion = "stars",
  lambda.log = TRUE,
  lambda.min.ratio = 1e-2,
  nlambda = 20,
  pulsar.select = TRUE,
  pulsar.params = list(
    thresh = 0.05,
    subsample.ratio = 0.8,
    rep.num = 50,
    seed = 1234
  ) )


getOptInd(t1)
opt <- getOptInd(t1)
nlam <- length(t1$lambda)
cat("Optimal lambda index:", opt, "\n")
cat("Number of lambda values:", nlam, "\n")
cat("Optimal lambda:", t1$lambda[opt], "\n")
cat("Minimum lambda:", min(t1$lambda), "\n")

if (opt == nlam) {
  cat("WARNING: Optimal lambda is at the minimum end of the lambda path.\n")
} else {
  cat("Optimal lambda is within the lambda path.\n")
}


adj <- getRefit(t1)

cat(
  "Number of nodes:",
  nrow(adj),
  "\n"
)

cat(
  "Number of edges:",
  Matrix::nnzero(adj) / 2,
  "\n"
)

cat(
  "Network sparsity:",
  sum(adj[upper.tri(adj)] != 0) /
    (nrow(adj) * (nrow(adj) - 1) / 2),
  "\n" )

getOptLambda(t1)
getOptInd(t1)

node_names <- taxa_names(data_f)

g <- graph_from_adjacency_matrix(
  adj, mode = "undirected", diag = FALSE )

# Assign OTU names
V(g)$name <- node_names


g <- delete_vertices(
  g,
  which(igraph::degree(g) == 0)
)

tax_sub <- tax_f[V(g)$name, , drop = FALSE]

# Check that all node names match taxonomy
stopifnot(
  all(V(g)$name %in% rownames(tax_f))
)

V(g)$Kingdom <- tax_sub[V(g)$name, "Kingdom"]
V(g)$Phylum  <- tax_sub[V(g)$name, "Phylum"]
V(g)$Class   <- tax_sub[V(g)$name, "Class"]
V(g)$Order   <- tax_sub[V(g)$name, "Order"]
V(g)$Family  <- tax_sub[V(g)$name, "Family"]
V(g)$Genus   <- tax_sub[V(g)$name, "Genus"]
V(g)$Species <- tax_sub[V(g)$name, "Species"]

adjacency_otu <- as.matrix(
  as_adjacency_matrix(
    g,
    sparse = FALSE
  )
)

write.csv(
  adjacency_otu,
  "1_0917_adjacency_raw_118_spie.csv",
  row.names = TRUE )

write_graph( g, file = "1_0917_raw_network_118_spie.graphml", format = "graphml" )

cat(
  "Final network nodes:",
  vcount(g),
  "\n"
)

cat(
  "Final network edges:",
  ecount(g),
  "\n"
)

cat(
  "Network density:",
  edge_density(g),
  "\n"
)


########## 模块提取 
otu_a = read.csv("field_data_386otu_ra_118.csv", header=TRUE, row.names=1) 
rownames(otu_a) <- otu_a[ , 12 ]
otu = otu_a[ , -c(1:12)]
otu_t = t(otu)
head(otu_t)
module0 = read.csv("1_0917_spie_raw_178_501_module.csv", header=TRUE, row.names=1)

# Extract the OTU names from module0
otu_list <- module0$v_name

# Extract the corresponding OTU row from otu_t
otu_module <- otu_t[rownames(otu_t) %in% otu_list, ]

write.csv( otu_module , "1_0917_spie_raw_178_501_module0.csv", row.names = T)

####
module = read.csv("1_0917_spie_raw_178_501_module1.csv", header=TRUE, row.names=1)

module_by_species <- aggregate(. ~ modularity_class, data = module , FUN = sum , na.rm = TRUE)

head(module_by_species[, 1:5])

write.csv(module_by_species, "1_0917_spie_raw_178_501_module2.csv", row.names = FALSE)



######## CA --mean cal-- 118 28species ----------------
library(plyr) 
library(ggtree)
library(ape) 
library(picante) 
library(phylosignal) 
library(phylobase) 
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_28sp")
data = read.csv("field_data_118.csv",header=TRUE, row.names=1)

data_mean <- ddply(data, c("Plant.species","Plant.species0",  "Plant.family"),summarise,   
                   P.concentrations = mean(P.concentrations,na.rm=TRUE),
                   Root.diameter = mean(Root.diameter,na.rm=TRUE),
                   Root.length.colonization = mean(Root.length.colonization,na.rm=TRUE),
                   
                   Glomeraceae = mean(Glomeraceae,na.rm=TRUE),
                   Acaulosporaceae = mean(Acaulosporaceae,na.rm=TRUE),
                   
                   Dominikia = mean(Dominikia, na.rm=TRUE),
                   Rhizophagus = mean(Rhizophagus,na.rm=TRUE),
                   
                   
                   Family.richness = mean(Family.richness,na.rm=TRUE),
                   Genus.richness = mean(Genus.richness,na.rm=TRUE),
                   OTU.richness = mean(OTU.richness, na.rm=TRUE),
                   
                   Module.1 = mean(Module.1, na.rm=TRUE),
                   Module.2 = mean(Module.2, na.rm=TRUE),
                   Module.3 = mean( Module.3 , na.rm=TRUE),
                   Module.4 = mean( Module.4 , na.rm=TRUE),
                   Module.5 = mean( Module.5 , na.rm=TRUE)
                   
)

data_mean

data.table::fwrite(data_mean,"field_data_118_28sp_mean.csv",  row.names = TRUE, bom = TRUE)


######## CA  v.phylomaker tree  28sp  --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_28sp")
df <- read.csv("28sp.csv")
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


######## CA -- pcoa rare 28 sp ------------------------------------
library(vegan)
library(ggplot2)
library(ggrepel)
library(ape)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_28sp")
data = read.csv("field_data_118_28sp.csv",header=TRUE, row.names=1) 
otu_a <- read.csv("field_data_386otu_rare_118_28sp.csv",header=TRUE, row.names=1)
otu = otu_a[, -c(1,2,3)] 
rownames(otu) = otu_a[, 2]
###
#otu_hellinger <- decostand(otu, method = "hellinger")

vare.dist <- vegdist(otu , method="bray")

pcoa_res <- pcoa(vare.dist)

pcoa_axis <- as.data.frame(pcoa_res$vectors[, 1:2])

colnames(pcoa_axis) <- c("PCoA1", "PCoA2")

eig_percent <- round(pcoa_res$values$Relative_eig[1:2] * 100, 2)
eig_percent

write.csv( pcoa_axis , '1_1025_rare_118_28sp_pcoa12.csv')  

###
set.seed(123)
# envfit
fit <- envfit( pcoa_res$vectors , data$Root.diameter, permutations = 9999)
fit
fit <- envfit( pcoa_res$vectors , data$P.concentrations, permutations = 9999)
fit

adonis2(vare.dist ~ data$Root.diameter, permutations = 9999) 
adonis2(vare.dist ~ data$P.concentrations, permutations = 9999) 
#  Mantel test
rd_dist <- dist(data$P.concentrations, method = "euclidean")
set.seed(123)
mantel_rd <- mantel(  vare.dist,  rd_dist,  method = "pearson",  permutations = 9999 )
mantel_rd

fit <- envfit( pcoa_res$vectors, data[, c("Root.diameter", "P.concentrations")], permutations = 9999 ) 

# 
fit_df <- data.frame( Axis1 = fit$vectors$arrows[, 1], Axis2 = fit$vectors$arrows[, 2],
                      r2 = fit$vectors$r,
                      P = fit$vectors$pvals )

fit_df$Variable <- rownames(fit$vectors$arrows)

# 调整箭头长度（可根据图形再微调）
arrow_mul <- 0.6
fit_df$PCoA1 <- fit_df$Axis1 * arrow_mul
fit_df$PCoA2 <- fit_df$Axis2 * arrow_mul
fit_df$Variable <- c("Root diameter", "Tissue P concentration")

# 准备数据
pcoa_axis$Sample <- rownames(otu_a)
eig_percent_formatted <- sprintf("%.2f", eig_percent)

# 基础 PCoA 图
pcoa <- ggplot(pcoa_axis, aes(x = PCoA1, y = PCoA2, )) + 
  geom_point( size = 2.5,  shape = 16,  fill = "grey66" ,color = "grey66", show.legend = F ) + scale_size(range = c(1, 6)) +
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


######## CA -- Table S4 -- Relationship between RD and AM fungal variables ----------------------
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
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_28sp")
data = read.csv("field_data_118_28sp.csv",header=TRUE, row.names=1) 
colnames(data)

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

figs7a <- ggplot( data , aes(y = Rhizophagus , x = Root.diameter  ))+   labs(title = "Field: 28 species") +    
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = expression("RA of "*italic(Rhizophagus)*" (%)")) +
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22",  show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth= 0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci1,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model11 ), linewidth = 0.8 , color = "#CB4641", linetype = "solid" ) + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "grey22", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 1.3, color = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.375, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.03, vjust = 2.5, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.375, 0.001 ), size = 4, inherit.aes = FALSE )
figs7a


##
lm_model <- lm( OTU.richness ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  OTU.richness ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model2 <- predict(pgls_model, newdata = data )  

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
pgls_ci2 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci2)

figs7b <- ggplot( data , aes(y = OTU.richness , x = Root.diameter  ))+      labs(title = "Field: 28 species") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "OTU richness" )+ 
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22",  show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid"  )+ 
  #geom_ribbon(data = pgls_ci2 ,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model2 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = - 4.5, colour = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.013, 0.560 ), size = 4 , inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = - 3, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.013, 0.560 ), size = 4 , inherit.aes = FALSE )
figs7b


##
lm_model <- lm( Family.richness ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  Family.richness ~ Root.diameter , data = comp_data, lambda = "ML") 
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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( Family.richness ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci3 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci3)

figs7c <- ggplot( data , aes(y = Family.richness , x = Root.diameter  ))+      labs(title = "Field: 28 species") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "Family richness" ) +
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22",  show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid"  )+ 
  #geom_ribbon(data = pgls_ci3,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model3 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = 1.1 , colour = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.005, 0.713 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = 2.2 , color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.005, 0.713 ), size = 4, inherit.aes = FALSE )
figs7c


##
lm_model <- lm( PCoA1 ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  PCoA1 ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model4 <- predict(pgls_model, newdata = data )  

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
pgls_ci4 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci4)

figs7d <- ggplot( data , aes( y = PCoA1 , x = Root.diameter  ))+      labs(title = "Field: 28 species") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "PCoA axis1 of AM fungal community" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22",  show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "black",  fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "dashed"  )+ 
  #geom_ribbon(data = pgls_ci4,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model4 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.1, vjust =  8.5 , colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.008, 0.655 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.1, vjust =  9.4, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.008, 0.655 ), size = 4, inherit.aes = FALSE )
figs7d


##
lm_model <- lm( PCoA2 ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  PCoA2 ~ Root.diameter , data = comp_data, lambda = "ML") 
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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( PCoA2 ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci5 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci5)

figs7e <- ggplot( data , aes( y = PCoA2 , x = Root.diameter  ))+      labs(title = "Field: 28 species") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "PCoA axis2 of AM fungal community" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22",  show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black",  fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci5,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model5 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.1, vjust =  2 , colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.224, 0.011 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.1, vjust =  3.2, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.130, 0.059 ), size = 4, inherit.aes = FALSE )
figs7e


##
lm_model <- lm( Module.3 ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  Module.3 ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model6 <- predict(pgls_model, newdata = data )  

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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( Module.3 ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci6 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci6)

figs7f <- ggplot( data , aes( y = Module.3 , x = Root.diameter ))+       labs(title = "Field: 28 species") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "RA of Module #3 (%)" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci6,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model6 ), linewidth = 0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.06, vjust = 1.3, colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.284, 0.004 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 2.5, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.284, 0.004 ), size = 4, inherit.aes = FALSE )
figs7f


#
figs7 <- cowplot::plot_grid( figs7a,  figs7b, figs7c, 
                             figs7d, figs7e, figs7f,
                             
                             ncol = 3, nrow = 2, labels = c( "(a)", "(b)",'(c)',
                                                             "(d)", "(e)","(f)",  
                                                             "(g)", "(h)", "(i)",
                                                             "(j)", "(k)", "(l)" ), 
                             label_fontface = "bold", label_x = -0.02, label_y = 1.01 , align = c("hv"))
figs7

setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom")
#eoffice::topptx( figs7 , width = 9 ,  height = 7 , "2_figs7_118_28sp_9*7.pptx" ) 


####
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


figs8a <- ggplot( data, aes(y = P.concentrations , x = Rhizophagus  ))+       labs(title = "Field: 28 species") +
  scale_x_continuous(name = expression("RA of "*italic(Rhizophagus)*" (%)")) + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci2,  aes(x = Rhizophagus, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model22 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  #guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 1.3, colour = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.169 , 0.030 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 2.5, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.117 , 0.074 ), size = 4, inherit.aes = FALSE )

figs8a


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

figs8b <- ggplot( data, aes(y = P.concentrations , x = OTU.richness  ))+     labs(title = "Field: 28 species") +  
  scale_x_continuous(name ="OTU richness") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22",  show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid" )+ 
  #geom_ribbon(data = pgls_ci4,  aes(x = OTU.richness, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model4 ), linewidth = 0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33",  face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 4.2, colour = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'<'~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.001, 0.972 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 5.4, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.002, 0.843 ), size = 4, inherit.aes = FALSE )

figs8b


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


figs8c <- ggplot( data, aes(y = P.concentrations ,  x = Family.richness  )) +    labs(title = "Field: 28 species") +   
  scale_x_continuous(name ="Family richness") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci6,  aes(x = Family.richness, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model6 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5,color= "grey33",  face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust =  3.7, colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.161, 0.034 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.03, vjust =  4.7, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.154, 0.039 ), size = 4, inherit.aes = FALSE )

figs8c


###
lm_model <- lm( P.concentrations ~ PCoA1, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ PCoA1 , data = comp_data, lambda = "ML") 
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
pred_grid <- seq(min(dat0$PCoA1), max(dat0$PCoA1), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( P.concentrations ~ PCoA1 ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( PCoA1 = pred_grid) ) } } }

# CI
pgls_ci8 <- data.frame( PCoA1 = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


figs8d <- ggplot( data, aes( y = P.concentrations , x = PCoA1  )) +       labs(title = "Field: 28 species") +
  scale_x_continuous(name ="PCoA axis1 of AM fungal community") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22",  show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black",  fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "dashed" )+ 
  geom_ribbon(data = pgls_ci8,  aes(x = PCoA1, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model8 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text = element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust =  3.8, colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.138, 0.052 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust =  4.8, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.280, 0.004 ), size = 4, inherit.aes = FALSE )


figs8d


###
lm_model <- lm( P.concentrations ~ PCoA2, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ PCoA2 , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model9 <- predict(pgls_model, newdata = data )  

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
pgls_ci9 <- data.frame( PCoA2 = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


figs8e <- ggplot( data, aes(y = P.concentrations , x = PCoA2  )) +     labs(title = "Field: 28 species") +  
  scale_x_continuous(name ="PCoA axis2 of AM fungal community") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22",  show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black",  fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci9,  aes(x = PCoA2, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model9 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33",  face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text = element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust =  1.8, colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.315, 0.002 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust =  2.8, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.315, 0.002 ), size = 4, inherit.aes = FALSE )


figs8e


###
lm_model <- lm( P.concentrations ~ Module.3, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ Module.3 , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model222 <- predict(pgls_model, newdata = data )  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$Module.3), max(dat0$Module.3), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( P.concentrations ~ Module.3 ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Module.3 = pred_grid) ) } } }

# CI
pgls_ci10 <- data.frame( Module.3 = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


figs8f <- ggplot( data, aes(y = P.concentrations , x = Module.3  )) +   labs(title = "Field: 28 species") +    
  scale_x_continuous(name ="RA of Module #3 (%)") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22",  show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci10,  aes(x = Module.3, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model222 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 1.2, colour = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.360, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 2.3, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.327, 0.001 ), size = 4, inherit.aes = FALSE )

figs8f


figs8 <- cowplot::plot_grid( figs8a,  figs8b, figs8c, 
                             figs8d, figs8e, figs8f,
                             
                             ncol = 3, nrow = 2, labels = c( "(a)", "(b)",'(c)',
                                                             "(d)", "(e)","(f)",  
                                                             "(g)", "(h)", "(i)",
                                                             "(j)", "(k)", "(l)" ), 
                             label_fontface = "bold", label_x = -0.02, label_y = 1.01 , align = c("hv"))
figs8


setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom")
#eoffice::topptx( figs8 , width = 9 ,  height = 7 , "2_figs8_118_28sp_9*7.pptx" ) 


#########################################################################################################
######### Sensitivity analysis2 29 sp Sensitivity analysis between the all-species dataset and the typical-AM-species dataset ##################################################################
######## D3 Reads table mean  & Alpha_diversity & RA calculation ---------------
library(amplicon)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_29sp")
otu_all <- read.csv('field_data_389otu_raw_128.csv', header=T,row.names = 1)
colnames(otu_all)

# OTU data merged by species (average calculated)
otu2 = otu_all[ , -c(1:4, 7:12 )]

otu_by_species <- aggregate(. ~ Plant.species, 
                            data = otu2[, !names(otu2) %in% "Plant.family"], 
                            FUN = mean, na.rm = TRUE)

family_info <- unique(otu2[, c("Plant.species", "Plant.family")])
otu_by_species <- merge(otu_by_species, family_info, by = "Plant.species", all.x = TRUE)


otu_by_species <- otu_by_species[, c("Plant.species", "Plant.family", 
                                     setdiff(names(otu_by_species), 
                                             c("Plant.species", "Plant.family")))]
head(otu_by_species[, 1:5])

data.table::fwrite( otu_by_species , "field_data_389otu_raw_128_29sp.csv",  row.names = TRUE, bom = TRUE)

### 
otu = otu_by_species
rownames(otu) <- otu[, 1]
otu <- otu[ , -c(1,2)]

otu_int <- round(otu)
##2_rarecurve_6_6
rarecurve(otu_int,  step = 10,  cex.lab = 1.2,  main = NULL,
          label = FALSE,  lines = TRUE, col = "tomato", lwd = 0.8,
          xlab = "Number of reads per plant species",  ylab = "Number of OTU" )


#### RA calculation
Tax <- read.csv('field_data_389otu_taxon.csv', header=T, row.names = 1)

####
otu_all <- read.csv('field_data_389otu_raw_128.csv', header=T,row.names = 1)
colnames(otu_all)
rownames(otu_all) <- otu_all[, 12]
otu_t <- t(otu_all[ , -c(1:12)])

format2stamp(otu_t, Tax, thre = 0.0, prefix = "RA_128")



######## D3 Co-occurrence network analysis SPIEC-EASI -- 128samples raw ---------------------------
#install.packages("glmnet")
#setwd("/Users/zengjing/Downloads") 
#install.packages("SpiecEasi-1.1.1.tar.gz", repos = NULL, type = "source")
# Load required libraries
library(SpiecEasi)
library(phyloseq)
library(igraph)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_29sp")
# Load OTU, taxonomic, and sample data
otu_all <- read.csv("field_data_389otu_raw_128.csv",header=TRUE, row.names=1) 
rownames(otu_all) <- otu_all[, 12]
otu <- as.matrix( t(otu_all[,-c(1:12)]) )
tax <- read.csv("field_data_389otu_taxon.csv", header = T, row.names = 1)
sample <- otu_all
rownames(sample) <- sample[, 12]

prev <- rowSums(otu > 0) / ncol(otu)
summary(prev)

table(prev >= 0.05)
table(prev >= 0.10)
table(prev >= 0.15)
table(prev >= 0.20)

c( "all" = sum(prev > 0),
   "5%" = sum(prev >= 0.05),
   "10%" = sum(prev >= 0.10),
   "15%" = sum(prev >= 0.15),
   "20%" = sum(prev >= 0.20)  )

# Create phyloseq object
otu_f <- otu[prev >= 0.05, ]
dim(otu_f)

tax_f <- tax[rownames(otu_f), , drop = FALSE]
all(rownames(otu_f) == rownames(tax_f))

data_f <- phyloseq(
  otu_table(
    otu_f,
    taxa_are_rows = TRUE
  ),
  tax_table(
    as.matrix(tax_f)
  ),
  sample_data(sample)
)

data_f

set.seed(1234)

t1 <- spiec.easi(
  data_f,
  method = "glasso",
  sel.criterion = "stars",
  lambda.log = TRUE,
  lambda.min.ratio = 1e-2,
  nlambda = 20,
  pulsar.select = TRUE,
  pulsar.params = list(
    thresh = 0.05,
    subsample.ratio = 0.8,
    rep.num = 50,
    seed = 1234
  ) )

getOptInd(t1)
opt <- getOptInd(t1)
nlam <- length(t1$lambda)
cat("Optimal lambda index:", opt, "\n")
cat("Number of lambda values:", nlam, "\n")
cat("Optimal lambda:", t1$lambda[opt], "\n")
cat("Minimum lambda:", min(t1$lambda), "\n")

if (opt == nlam) {
  cat("WARNING: Optimal lambda is at the minimum end of the lambda path.\n")
} else {
  cat("Optimal lambda is within the lambda path.\n")
}


adj <- getRefit(t1)

cat(
  "Number of nodes:",
  nrow(adj),
  "\n"
)

cat(
  "Number of edges:",
  Matrix::nnzero(adj) / 2,
  "\n"
)

cat(
  "Network sparsity:",
  sum(adj[upper.tri(adj)] != 0) /
    (nrow(adj) * (nrow(adj) - 1) / 2),
  "\n" )

getOptLambda(t1)
getOptInd(t1)

node_names <- taxa_names(data_f)

g <- graph_from_adjacency_matrix(
  adj, mode = "undirected", diag = FALSE )

# Assign OTU names
V(g)$name <- node_names


g <- delete_vertices(
  g,
  which(igraph::degree(g) == 0)
)

tax_sub <- tax_f[V(g)$name, , drop = FALSE]

# Check that all node names match taxonomy
stopifnot(
  all(V(g)$name %in% rownames(tax_f))
)

V(g)$Kingdom <- tax_sub[V(g)$name, "Kingdom"]
V(g)$Phylum  <- tax_sub[V(g)$name, "Phylum"]
V(g)$Class   <- tax_sub[V(g)$name, "Class"]
V(g)$Order   <- tax_sub[V(g)$name, "Order"]
V(g)$Family  <- tax_sub[V(g)$name, "Family"]
V(g)$Genus   <- tax_sub[V(g)$name, "Genus"]
V(g)$Species <- tax_sub[V(g)$name, "Species"]

adjacency_otu <- as.matrix(
  as_adjacency_matrix(
    g,
    sparse = FALSE
  )
)

write.csv(
  adjacency_otu,
  "1_1024_adjacency_raw_128_spie.csv",
  row.names = TRUE )

write_graph( g, file = "1_1024_raw_network_128_spie.graphml", format = "graphml" )

cat(
  "Final network nodes:",
  vcount(g),
  "\n"
)

cat(
  "Final network edges:",
  ecount(g),
  "\n"
)

cat(
  "Network density:",
  edge_density(g),
  "\n"
)



########## 模块提取 
otu_a = read.csv("field_data_389otu_ra_128.csv", header=TRUE, row.names=1) 
rownames(otu_a) <- otu_a[ , 12 ]
otu = otu_a[ , -c(1:12)]
otu_t = t(otu)
head(otu_t)
module0 = read.csv("1_1024_spie_raw_129_239_module.csv", header=TRUE, row.names=1)

# Extract the OTU names from module0
otu_list <- module0$v_name

# Extract the corresponding OTU row from otu_t
otu_module <- otu_t[rownames(otu_t) %in% otu_list, ]

write.csv( otu_module , "1_1024_spie_raw_129_239_module0.csv", row.names = T)

####
module = read.csv("1_1024_spie_raw_129_239_module1.csv", header=TRUE, row.names=1)

module_by_species <- aggregate(. ~ modularity_class, data = module , FUN = sum , na.rm = TRUE)

head(module_by_species[, 1:5])

write.csv(module_by_species, "1_1024_spie_raw_129_239_module2.csv", row.names = FALSE)



######## D3 --mean cal -- 128 29species ----------------
library(plyr)
library(ggtree)
library(ape) 
library(picante) 
library(phylosignal) 
library(phylobase) 
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_29sp")
data = read.csv("field_data_128.csv",header=TRUE, row.names=1)

data_mean <- ddply(data, c("Plant.species","Plant.species0",  "Plant.family"),summarise,   
                   P.concentrations = mean(P.concentrations,na.rm=TRUE),
                   Root.diameter = mean(Root.diameter,na.rm=TRUE),
                   Root.length.colonization = mean(Root.length.colonization,na.rm=TRUE),
                   
                   Glomeraceae = mean(Glomeraceae,na.rm=TRUE),
                   Acaulosporaceae = mean(Acaulosporaceae,na.rm=TRUE),
                   
                   Dominikia = mean(Dominikia, na.rm=TRUE),
                   Rhizophagus = mean(Rhizophagus,na.rm=TRUE),
                   
                   
                   Family.richness = mean(Family.richness,na.rm=TRUE),
                   Genus.richness = mean(Genus.richness,na.rm=TRUE),
                   OTU.richness = mean(OTU.richness, na.rm=TRUE),
                   
                   Module.1 = mean(Module.1, na.rm=TRUE),
                   Module.2 = mean(Module.2, na.rm=TRUE),
                   Module.3 = mean( Module.3 , na.rm=TRUE),
                   Module.4 = mean( Module.4 , na.rm=TRUE),
                   Module.5 = mean( Module.5 , na.rm=TRUE)
                   
)

data_mean

data.table::fwrite(data_mean,"field_data_128_29sp_mean.csv",  row.names = TRUE, bom = TRUE)


######## D3 v.phylomaker tree  29 sp  --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_29sp")
df <- read.csv("29sp.csv")
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

write.tree(sp_tr, "1_29sp_phylomakers1.nwk")

#sp_tr$node.label <- NULL
#write.tree(sp_tr, "1-35sp-phylomakers1-mega.nwk")

is.rooted(sp_tr)       
is.binary.tree(sp_tr)  

####
library(dplyr)
tree <- read.tree("1_29sp_phylomakers1.nwk")
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

write.tree( tree, file = "1_29sp_phylomakers1_abbr.nwk")


######## D3 --  -- pcoa 29 sp ------------------------------------
library(vegan)
library(ggplot2)
library(ggrepel)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_29sp")
data = read.csv("field_data_128_29sp.csv",header=TRUE, row.names=1) 
otu_a <- read.csv("field_data_389otu_raw_128_29sp.csv",header=TRUE, row.names=1)
otu = otu_a[, -c(1,2)] 
rownames(otu) = otu_a[, 1]
###
#otu_hellinger <- decostand(otu, method = "hellinger")

vare.dist <- vegdist(otu , method="bray")

pcoa_res <- pcoa(vare.dist)

pcoa_axis <- as.data.frame(pcoa_res$vectors[, 1:2])

colnames(pcoa_axis) <- c("PCoA1", "PCoA2")

eig_percent <- round(pcoa_res$values$Relative_eig[1:2] * 100, 2)
eig_percent

write.csv( pcoa_axis , '1_1024_rare_128_29sp_pcoa12.csv')  

###
set.seed(123)
# envfit
fit <- envfit( pcoa_res$vectors , data$Root.diameter, permutations = 9999)
fit
fit <- envfit( pcoa_res$vectors , data$P.concentrations, permutations = 9999)
fit

adonis2(vare.dist ~ data$Root.diameter, permutations = 9999) 
adonis2(vare.dist ~ data$P.concentrations, permutations = 9999) 
#  Mantel test
rd_dist <- dist(data$P.concentrations, method = "euclidean")
set.seed(123)
mantel_rd <- mantel(  vare.dist,  rd_dist,  method = "pearson",  permutations = 9999 )
mantel_rd

fit <- envfit( pcoa_res$vectors, data[, c("Root.diameter", "P.concentrations")], permutations = 9999 ) 

# 
fit_df <- data.frame( Axis1 = fit$vectors$arrows[, 1], Axis2 = fit$vectors$arrows[, 2],
                      r2 = fit$vectors$r,
                      P = fit$vectors$pvals )

fit_df$Variable <- rownames(fit$vectors$arrows)

# 调整箭头长度（可根据图形再微调）
arrow_mul <- 0.6
fit_df$PCoA1 <- fit_df$Axis1 * arrow_mul
fit_df$PCoA2 <- fit_df$Axis2 * arrow_mul
fit_df$Variable <- c("Root diameter", "Tissue P concentration")

# 准备数据
pcoa_axis$Sample <- rownames(otu_a)
eig_percent_formatted <- sprintf("%.2f", eig_percent)

# 基础 PCoA 图
pcoa <- ggplot(pcoa_axis, aes(x = PCoA1, y = PCoA2, )) + 
  geom_point( size = 2.5,  shape = 16,  fill = "grey66" ,color = "grey66", show.legend = F ) + scale_size(range = c(1, 6)) +
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


######## D3 -- Table S3 -- Relationship between RD and AM fungal variables ----------------------
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
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/2_amfcom_29sp")
data = read.csv("field_data_128_29sp.csv",header=TRUE, row.names=1) 
colnames(data)

##
rooted_tree <- read.tree("1_29sp_phylomakers1_abbr.nwk") 
plot(rooted_tree)

rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data, names.col = "Plant.species", vcv = T, vcv.dim = 2 )
nrow(comp_data$data)
nrow(data)
setdiff(data$Plant.species, rownames(comp_data$data))

##
lm_model <- lm( Rhizophagus ~ Root.diameter, data = data ) 
summary(lm_model)
pgls_model <- pgls(  Rhizophagus ~ Root.diameter, data = comp_data, lambda = "ML" )  
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

fig2a <- ggplot( data , aes(y = Rhizophagus , x = Root.diameter  ))+   labs(title = "Field: 29 species") +    
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = expression("RA of "*italic(Rhizophagus)*" (%)")) +
  geom_point( size = 2.5,  shape = 16,  fill = "grey66" ,color = "grey66", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth= 0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci1,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model11 ), linewidth = 0.8 , color = "#CB4641", linetype = "solid" ) + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "grey22", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 1.3, color = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.375, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.03, vjust = 2.5, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.375, 0.001 ), size = 4, inherit.aes = FALSE )
fig2a


##
lm_model <- lm( OTU.richness ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  OTU.richness ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model2 <- predict(pgls_model, newdata = data )  

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
pgls_ci2 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci2)

fig2b <- ggplot( data , aes(y = OTU.richness , x = Root.diameter  ))+      labs(title = "Field: 29 species") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "OTU richness" )+ 
  geom_point( size = 2.5,  shape = 16,  fill = "grey66" ,color = "grey66", show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid"  )+ 
  #geom_ribbon(data = pgls_ci2 ,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model2 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = - 4.5, colour = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.013, 0.560 ), size = 4 , inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = - 3, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.013, 0.560 ), size = 4 , inherit.aes = FALSE )
fig2b


##
lm_model <- lm( Family.richness ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  Family.richness ~ Root.diameter , data = comp_data, lambda = "ML") 
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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( Family.richness ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci3 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci3)

fig2c <- ggplot( data , aes(y = Family.richness , x = Root.diameter  ))+      labs(title = "Field: 29 species") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "Family richness" ) +
  geom_point( size = 2.5,  shape = 16,  fill = "grey66" ,color = "grey66", show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid"  )+ 
  #geom_ribbon(data = pgls_ci3,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model3 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = 1.1 , colour = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.005, 0.713 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = 2.2 , color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.005, 0.713 ), size = 4, inherit.aes = FALSE )
fig2c


##
lm_model <- lm( PCoA1 ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  PCoA1 ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model4 <- predict(pgls_model, newdata = data )  

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
pgls_ci4 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci4)

fig2d <- ggplot( data , aes( y = PCoA1 , x = Root.diameter  ))+      labs(title = "Field: 29 species") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "PCoA axis1 of AM fungal community" )+
  geom_point( size = 2.5,  shape = 16,  fill = "grey66" ,color = "grey66", show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "black",  fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "dashed"  )+ 
  #geom_ribbon(data = pgls_ci4,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model4 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.1, vjust =  8.5 , colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.008, 0.655 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.1, vjust =  9.4, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.008, 0.655 ), size = 4, inherit.aes = FALSE )
fig2d


##
lm_model <- lm( PCoA2 ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  PCoA2 ~ Root.diameter , data = comp_data, lambda = "ML") 
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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( PCoA2 ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci5 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci5)

fig2e <- ggplot( data , aes( y = PCoA2 , x = Root.diameter  ))+      labs(title = "Field: 29 species") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "PCoA axis2 of AM fungal community" )+
  geom_point( size = 2.5,  shape = 16,  fill = "grey66" ,color = "grey66",  show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black",  fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci5,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model5 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.1, vjust =  2 , colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.224, 0.011 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.1, vjust =  3.2, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.130, 0.059 ), size = 4, inherit.aes = FALSE )
fig2e


##
lm_model <- lm( Module.2 ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  Module.2 ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model6 <- predict(pgls_model, newdata = data )  

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
pgls_ci6 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci6)

fig2f <- ggplot( data , aes( y = Module.2 , x = Root.diameter ))+       labs(title = "Field: 29 species") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "RA of Module #3 (%)" )+
  geom_point( size = 2.5,  shape = 16,  fill = "grey66" ,color = "grey66",  show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci6,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model6 ), linewidth = 0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.06, vjust = 1.3, colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.284, 0.004 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 2.5, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.284, 0.004 ), size = 4, inherit.aes = FALSE )
fig2f



####
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


fig3a <- ggplot( data, aes(y = P.concentrations , x = Rhizophagus  ))+       labs(title = "Field: 29 species") +
  scale_x_continuous(name = expression("RA of "*italic(Rhizophagus)*" (%)")) + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 16,  fill = "grey66" ,color = "grey66", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci2,  aes(x = Rhizophagus, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model22 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  #guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 1.3, colour = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.169 , 0.030 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 2.5, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.117 , 0.074 ), size = 4, inherit.aes = FALSE )

fig3a


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

fig3b <- ggplot( data, aes(y = P.concentrations , x = OTU.richness  ))+     labs(title = "Field: 29 species") +  
  scale_x_continuous(name ="OTU richness") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 16, fill = "grey66" ,color = "grey66", show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid" )+ 
  #geom_ribbon(data = pgls_ci4,  aes(x = OTU.richness, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model4 ), linewidth = 0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33",  face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 4.2, colour = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'<'~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.001, 0.972 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 5.4, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.002, 0.843 ), size = 4, inherit.aes = FALSE )

fig3b


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


fig3c <- ggplot( data, aes(y = P.concentrations ,  x = Family.richness  )) +    labs(title = "Field: 29 species") +   
  scale_x_continuous(name ="Family richness") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 16,  fill = "grey66" ,color = "grey66", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci6,  aes(x = Family.richness, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model6 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5,color= "grey33",  face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust =  3.7, colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.161, 0.034 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.03, vjust =  4.7, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.154, 0.039 ), size = 4, inherit.aes = FALSE )

fig3c


###
lm_model <- lm( P.concentrations ~ PCoA1, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ PCoA1 , data = comp_data, lambda = "ML") 
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
pred_grid <- seq(min(dat0$PCoA1), max(dat0$PCoA1), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( P.concentrations ~ PCoA1 ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( PCoA1 = pred_grid) ) } } }

# CI
pgls_ci8 <- data.frame( PCoA1 = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig3d <- ggplot( data, aes( y = P.concentrations , x = PCoA1  )) +       labs(title = "Field: 29 species") +
  scale_x_continuous(name ="PCoA axis1 of AM fungal community") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 16,  fill = "grey66" ,color = "grey66", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black",  fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "dashed" )+ 
  geom_ribbon(data = pgls_ci8,  aes(x = PCoA1, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model8 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text = element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust =  3.8, colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.138, 0.052 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust =  4.8, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.280, 0.004 ), size = 4, inherit.aes = FALSE )


fig3d


###
lm_model <- lm( P.concentrations ~ PCoA2, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ PCoA2 , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model9 <- predict(pgls_model, newdata = data )  

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
pgls_ci9 <- data.frame( PCoA2 = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig3e <- ggplot( data, aes(y = P.concentrations , x = PCoA2  )) +     labs(title = "Field: 29 species") +  
  scale_x_continuous(name ="PCoA axis2 of AM fungal community") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 16, fill = "grey66" ,color = "grey66", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black",  fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci9,  aes(x = PCoA2, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model9 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33",  face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text = element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust =  1.8, colour = "black", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.315, 0.002 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust =  2.8, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.315, 0.002 ), size = 4, inherit.aes = FALSE )


fig3e


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


fig3f <- ggplot( data, aes(y = P.concentrations , x = Module.2  )) +   labs(title = "Field: 29 species") +    
  scale_x_continuous(name ="RA of Module #3 (%)") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 16,  fill = "grey66" ,color = "grey66", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci10,  aes(x = Module.2, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model222 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 1.2, colour = "black",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.360, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 2.3, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.327, 0.001 ), size = 4, inherit.aes = FALSE )

fig3f


#########################################################################################################
######### AM Fungal Community Samples – Bias Removal Test ---------------------
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom/amfcom_35sp")
data = read.csv("field_data_391otu_raw_155.csv", header=TRUE, row.names=1)
colnames(data)

library(ggplot2)
##
table(data$P20)

##
plot_data <- data[, c("Root.diameter", "P20")]
plot_data <- plot_data[complete.cases(plot_data), ]

# dele
dele_data <- plot_data[plot_data$P20 == "Dele", ]

### 
fig.r1 = ggplot(  plot_data, aes( x = Root.diameter,  y = P20 ) ) +
  geom_jitter( aes(color = P20), width = 0,  height = 0.12, size = 2.5,   alpha = 0.75 ) +
  scale_color_manual( values = c("Dele" = "#d95f02",  "Retained" = "grey70"    )  ) +
  labs(  x = "Root diameter (mm)",  y = NULL,  color = NULL ) + theme_classic() +
  theme(axis.text = element_text(size=12.5), axis.text.y = element_text(angle = 0, hjust = 0.5), 
        legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  legend.position = "none",
        axis.title=element_text(size=12.5)) 

fig.r1

setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/2_amfcom") 
eoffice::topptx( fig.r1 , width = 5 ,  height = 4 , "3_figr1.pptx" ) 

###
data$Deleted <- ifelse(data$P20 == "Dele", 1, 0)

model <- glm( Deleted ~ Root.diameter,  data = data,  family = binomial )
summary(model)

wilcox.test( Root.diameter ~ P20,  data = data)


###
data$Graminoid <- ifelse( data$Plant.family %in% c( "Poaceae", "Cyperaceae"), "Graminoid", "Non-graminoid" )
tab <- table(data$Graminoid, data$Deleted)
tab

fisher.test(tab)

table(data$Plant.family, data$Deleted)
tab <- table(data$Graminoid, data$Deleted)




