
######## OTU table rarefaction and diversity analysis ---------------
library(vegan)
library(picante)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
otu_all<- read.csv("field_data_otu_raw.csv",header=TRUE, row.names=1)

otu <- otu_all[,-c(1:4)]
otu_rrarefy = rrarefy(otu, median(rowSums(otu)))  

#write.csv( otu_rrarefy, "field_data_otu_rrarefy.csv")

Alpha_diversity_index <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Obs <-  est[1, ]
  Shannon <- vegan::diversity(x, index = 'shannon', base = base)
  Simpson <- vegan::diversity(x, index = 'simpson')
  Pielou <- Shannon / log(Obs, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  result <- rbind(est, Shannon, Simpson, Pielou, goods_coverage)
  
  if (!is.null(tree)) {
    Pd <- pd(x, tree, include.root = FALSE)[1]
    Pd <- t(Pd)
    result <- rbind(result, Pd)
  }
  
  result <- as.data.frame(t(result))
  return(result)
}

alpha_diversity <- Alpha_diversity_index(otu_rrarefy)

######## Co-occurrence network analysis ---------------
library(psych)
library(reshape2)
library(igraph)
library(Hmisc)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
otu = read.csv("field_data_otu_network.csv",header=TRUE, row.names=1) 

cor = rcorr(t(otu), type = 'spearman')  

r = cor$r 
p = cor$P                       
p = p.adjust(p, method = 'BH') 

r[p > 0.05 | abs(r) < 0.6] = 0  

g = graph_from_adjacency_matrix(r, mode ="undirected",weighted=TRUE,diag = FALSE)
g
g = delete_vertices(g, names(degree(g)[degree(g) == 0]))

E(g)$corr = E(g)$weight       
E(g)$weight = abs(E(g)$weight)   

tax = read.csv('field_data_otu_taxon.csv', row.names=1, header=T)    

tax = tax[as.character(V(g)$name), ]          
V(g)$Kingdom = tax$Kingdom                          
V(g)$Phylum = tax$Phylum                           
V(g)$Class = tax$Class                              
V(g)$Order = tax$Order                           
V(g)$Family = tax$Family                           
V(g)$Genus = tax$Genus                               
V(g)$Species = tax$Species                        

node_list = data.frame(
  label = names(V(g)),
  kingdom = V(g)$Kingdom,
  phylum = V(g)$Phylum,
  class = V(g)$Class,
  order = V(g)$Order,
  family = V(g)$Family,
  genus=V(g)$Genus,
  species = V(g)$Species )   

head(node_list)
#write.csv(node_list, '1_network_node.csv')   

edge = data.frame(as_edgelist(g)) 
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)

head(edge_list)
#write.csv(edge_list, '1_network_edge.csv')  

#write_graph(g, '1_field_data_otu_network.graphml', format = 'graphml')                     


######## Nestedness and modularity of plant-AM fungal bipartite network ---------------
library(bipartite) 
library(vegan)
library(dplyr) 
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
data = read.csv("field_data_otu_ra.csv", header=TRUE, row.names=1) 

sample_matrix <- t(data) 
head(sample_matrix, c(5, 5))
binary_matrix <- ifelse(sample_matrix > 0, 1, 0)

set.seed(1234)
null_model <- nullmodel( binary_matrix , N = 1000, method = "r2dtable")

nested_NODF <- function(x) nested(x, method = "NODF")

result <- oecosimu(
  binary_matrix,
  nested_NODF,
  method = "r2dtable",
  nsimul = 1000
)

result


####
mod_result <- computeModules(binary_matrix)

plotModuleWeb(mod_result)

mod_result@likelihood

mod_obs <- mod_result
Q_obs <- mod_obs@likelihood

n_sim <- 1000
Q_null <- numeric(n_sim)

for (i in 1:n_sim) {
  null_mat <- r2dtable(1, rowSums(binary_matrix), colSums(binary_matrix))[[1]]
  mod_null <- computeModules(null_mat)
  Q_null[i] <- mod_null@likelihood
}


p_value <- sum(Q_null >= Q_obs) / n_sim
p_value

hist(Q_null, main = "Null distribution of modularity", xlab = "Modularity (Q)")
abline(v = Q_obs, col = "red", lwd = 2)


######## Effects of plant species identity on plant traits and AM fungal variables ----------------
library(nlme)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
data = read.csv("field_data_all.csv",header=TRUE, row.names=1) 

fm <- lme(Root.diameter ~ Plant.species, data = data , random=~1|Block)
anova(fm)
fm <- lme(P.concentrations ~ Plant.species, data = data , random=~1|Block)
anova(fm)


fm <- lme(Root.length.colonization ~ Plant.species, data = data , random=~1|Block)
anova(fm)

fm <- lme(RA.of.Glomeraceae ~ Plant.species, data = data , random=~1|Block)
anova(fm)
fm <- lme(RA.of.Acaulosporaceae ~ Plant.species, data = data , random=~1|Block)
anova(fm)
fm <- lme(RA.of.Dominikia ~ Plant.species, data = data , random=~1|Block)
anova(fm)
fm <- lme(RA.of.Rhizophagus ~ Plant.species, data = data , random=~1|Block)
anova(fm)
fm <- lme(RA.of.OTU1 ~ Plant.species, data = data , random=~1|Block)
anova(fm)
fm <- lme(RA.of.OTU2 ~ Plant.species, data = data , random=~1|Block)
anova(fm)

fm <- lme(Family.Richness ~ Plant.species, data = data , random=~1|Block)
anova(fm)
fm <- lme(OTU.Richness ~ Plant.species, data = data , random=~1|Block)
anova(fm)

fm <- lme(RA.of.Module..1 ~ Plant.species, data = data , random=~1|Block)
anova(fm)


######## Plant phylogenetic signals of plant traits and AM fungal variables ----------------
library(plyr) 
library(ggtree)
library(ape)
library(picante)
library(phylosignal)
library(phylobase)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys")  
data = read.csv("field_data_all.csv",header=TRUE, row.names=1)

data_mean <- ddply(data, c("Plant.species"),summarise,   
                     P.concentrations = mean(P.concentrations,na.rm=TRUE),
                     Root.diameter = mean(Root.diameter,na.rm=TRUE),
                     Root.length.colonization = mean(Root.length.colonization,na.rm=TRUE),
                     Family.Richness = mean(Family.Richness, na.rm=TRUE), 
                     OTU.Richness = mean(OTU.Richness,na.rm=TRUE),
                     RA.of.Module..1 = mean(RA.of.Module..1,na.rm=TRUE),
                     RA.of.Glomeraceae = mean(RA.of.Glomeraceae,na.rm=TRUE),
                     RA.of.Acaulosporaceae = mean(RA.of.Acaulosporaceae,na.rm=TRUE),
                     RA.of.Dominikia = mean(RA.of.Dominikia,na.rm=TRUE),
                     RA.of.Rhizophagus = mean(RA.of.Rhizophagus,na.rm=TRUE),
                     RA.of.OTU1 = mean(RA.of.OTU1,na.rm=TRUE),
                     RA.of.OTU2 = mean(RA.of.OTU2,na.rm=TRUE)   )

data_mean

#write.csv(data_mean,"field_data_mean.csv")

data_mean1 = data_mean
rownames(data_mean1) <- data_mean[, 1]
data_mean1 <- data_mean1[, -1]

##
tree_raw <- read.tree("2021mqamf_35sp.nwk") 
plot(tree_raw)
rooted_tree1 <- root(tree_raw, outgroup = "Juncus", node = NULL, resolve.root = T, interactive = FALSE, edgelabel = FALSE)
plot(rooted_tree1)
rooted_tree  <- drop.tip(rooted_tree1, "Juncus")
plot(rooted_tree)


tree_p <- as(rooted_tree, "phylo4")
class(tree_p)
p4d <- phylo4d(tree_p, data_mean1)
phy.sig <- phyloSignal(p4d, reps = 9999)
phy.sig
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


cols <- c("Amaryllidaceae"= "#6C8FA6", "Cyperaceae" = "#5F97C6", "Poaceae" = "#8FB9D1", 
          "Ranunculaceae" = "#BFA87F","Onagraceae" = "#B86A2E", "Brassicaceae"= "#B09AA8",  
          "Rosaceae" =  "#E0A16F", "Fabaceae" = "#FFA64D" , "Caryophyllaceae" = "#7FC1B5" , 
          "Polygonaceae"= "#6BA89F",  "Gentianaceae"= "#9e8c93","Lamiaceae"  = "#D94F30",
          "Plantaginaceae" = "#A83248", "Apiaceae" = "#B85A78", "Asteraceae" = "#E9746A") 


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
  labs(x="",y="Tissue P concentrations (mg/g)")+
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


######## Relationship between plant traits and AM fungal colonization ---------------
library(ggplot2)
library(gridExtra)
library(reshape)
library(gam)
library(segmented)
library(chngpt)
library(ggpmisc)
library(car) 
library(doBy)
library(dplyr) 
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
data = read.csv("field_data_all.csv",header=TRUE, row.names=1) 

cols <- c("Amaryllidaceae"= "#6C8FA6", "Cyperaceae" = "#5F97C6", "Poaceae" = "#8FB9D1", 
          "Ranunculaceae" = "#BFA87F","Onagraceae" = "#B86A2E", "Brassicaceae"= "#B09AA8",  
          "Rosaceae" =  "#E0A16F", "Fabaceae" = "#FFA64D" , "Caryophyllaceae" = "#7FC1B5" , 
          "Polygonaceae"= "#6BA89F",  "Gentianaceae"= "#9e8c93","Lamiaceae"  = "#D94F30",
          "Plantaginaceae" = "#A83248", "Apiaceae" = "#B85A78", "Asteraceae" = "#E9746A") 

se <- function(x) sd(x)/sqrt(length(x))
conf_95 <- function(x) t.test(x, conf.level = 0.95)$conf.int[2]

outlier_test = outlierTest(lm( Root.length.colonization ~ Root.diameter , data = data)) 
outlier_test 
sumdata <- summaryBy( Root.diameter + Root.length.colonization ~ Plant.species, data, FUN = c(mean, sd, se, conf_95))
sumdata
sumdata$Family <- data$Plant.family[match(sumdata$Plant.species, data$Plant.species)]
sumdata$Root.diameter <- sumdata$Root.diameter.mean
sumdata$Root.length.colonization <- sumdata$Root.length.colonization.mean
colnames(sumdata)

fig1a <- ggplot(sumdata, aes(y = Root.length.colonization , x = Root.diameter, color = Family  ))+      
  scale_color_manual(values = cols ) +  labs(color="Plant family")+
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "Root length colonization (%)" )+
  geom_errorbarh(aes(xmin = Root.diameter + Root.diameter.se, xmax = Root.diameter - Root.diameter.se), height = 0, show.legend = F) + # color = "grey77", 
  geom_errorbar(aes(ymin = Root.length.colonization + Root.length.colonization.se, ymax = Root.length.colonization - Root.length.colonization.se), width = 0, show.legend = F) +
  geom_point(size =  2.5, alpha = 1, shape = 16, aes(color = Family), show.legend = F ) + 
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")), data = data , formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.98,label.y = 0.02, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  theme_bw() + theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=11), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))
fig1a


outlier_test = outlierTest(lm( P.concentrations ~ Root.length.colonization , data = data)) 
outlier_test 
sumdata <- summaryBy( P.concentrations + Root.length.colonization ~ Plant.species, data, FUN = c(mean, sd, se ))
sumdata
sumdata$Family <- data$Plant.family[match(sumdata$Plant.species, data$Plant.species)]
sumdata$P.concentrations <- sumdata$P.concentrations.mean
sumdata$Root.length.colonization <- sumdata$Root.length.colonization.mean
colnames(sumdata)

fig1b <- ggplot(sumdata, aes(y = P.concentrations , x = Root.length.colonization, color = Family  ))+      
  scale_x_continuous(name ="Root length colonization (%)") + scale_y_continuous(name = "TP (mg/g)" )+
  geom_errorbarh(aes(xmin = Root.length.colonization + Root.length.colonization.se, xmax = Root.length.colonization - Root.length.colonization.se), height = 0,show.legend = F) +
  geom_errorbar(aes(ymin = P.concentrations + P.concentrations.se, ymax = P.concentrations - P.concentrations.se), width = 0, show.legend = F ) +
  geom_point(size =  2.5, alpha = 1, shape = 16, aes(color = Family), show.legend = F ) + 
  scale_color_manual(values = cols, breaks = c("Asteraceae",  "Plantaginaceae", "Lamiaceae",
                                               "Gentianaceae",   "Rosaceae", "Caryophyllaceae",
                                               "Ranunculaceae", "Poaceae","Cyperaceae", "Amaryllidaceae" ) )+  labs(color="Plant family")+
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")), data = data, formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.02, label.y = 0.98, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  theme_bw() + theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=11), legend.title = element_text(size=12, hjust=0.15),  axis.title=element_text(size=12)) + #
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", label.hjust = 0))
fig1b 

#
legend <- cowplot::get_legend(fig1b) 

###
colnames(data)
dat <- data[ , c( 3, 5 , 6, 7)]
colnames(dat)
head(dat)

set.seed(12345)
bootari <- function(DF, model, var,
                    window = 85,
                    NsimB = 500,
                    ncore = max(1, parallel::detectCores() - 2)) {
  stopifnot(var %in% colnames(DF))
  xvar <- DF[[var]]
  if (is.data.frame(xvar)) xvar <- xvar[[1]]
  if (!is.numeric(xvar)) stop("var must be numeric")
  DF[[var]] <- xvar
  DF <- DF[order(DF[[var]]), ]
  vectorVAR <- sort(unique(DF[[var]]))
  nVAR <- length(vectorVAR)
  
  if (nVAR <= window) {
    stop("Not enough unique values in var for the chosen window size")
  }
  
  library(lme4)
  library(boot)
  library(doSNOW)
  
  cl <- makeCluster(ncore)
  registerDoSNOW(cl)
  
  bt <- NULL
  varNote <- numeric()
  
  for (i in 1:(nVAR - window)) {
    
    try({
      
      winVAR <- vectorVAR[i:(i + window)]
      
      dfi <- DF[DF[[var]] %in% winVAR, ]
      
      
      dfi$P.concentrations <-
        scale(dfi$P.concentrations, center = TRUE, scale = TRUE)
      dfi$Root.length.colonization <-
        scale(dfi$Root.length.colonization, center = TRUE, scale = TRUE)
      dfi$Root.diameter <- dfi$Root.diameter 
      
      fit <- lmer(as.formula(model), data = dfi, REML = FALSE)
      
      bt_i <- bootMer(
        fit,
        FUN = fixef,
        nsim = NsimB,
        parallel = "snow",
        cl = cl
      )$t
      
      
      bt <- rbind(bt, bt_i)
      varNote <- c(varNote, rep(mean(winVAR), NsimB))
      
    }, silent = TRUE)
  }
  
  stopCluster(cl)
  
  
  bt <- as.data.frame(bt)
  colnames(bt) <- names(fixef(fit))
  bt$var <- varNote
  
  return(bt)
}


model_1 <- "P.concentrations ~ Root.length.colonization + (1|Block)"
var <- "Root.diameter"

btaridity <- bootari(
  DF = dat,
  model = model_1,
  var = var,
  window = 85,
  NsimB = 500 )
colnames(btaridity) <- c("Intcp","TP","var")

model_11 <-"P.concentrations ~ Root.length.colonization+(1|Block)"
var1<-"Root.diameter"
btaridity1 <- bootari(
  DF = dat,
  model = model_11,
  var = var1,
  window = 85,
  NsimB = 500 )
colnames(btaridity1) <- c("Intcp","TP","var")

lsbt<-list(bt1=btaridity, bt2=btaridity1)
lsgood<-lsbt

for (n in 1:length(lsbt)){
  btn<-lsbt[[n]]
  btn<-btn[complete.cases(btn),]
  for (i in 1:length(unique(btn$var))){
    bti<-btn[btn$var==unique(btn$var)[i],]
    for(j in 1:(length(colnames(bti))-1)){
      x<-bti[,j]
      upci<-quantile(x, 0.95)
      lowci<-quantile(x,0.05)
      x[x>=upci | x<=lowci]<-NA
      bti[,j]<-x
    }
    btn[btn$var==unique(btn$var)[i],]<-bti
  }
  lsgood[[n]]<-btn
}

df1=lsbt[[1]]
colnames(df1)
length(unique(df1$var))
table(df1$var)

#write.csv(df1, "1_150_MW_85_TP_RLC_RAD.csv", quote=F)

df1=read.csv("1_150_MW_85_TP_RLC_RAD.csv", head=T, row.name = 1)

df_mean <- df1 %>%
  dplyr::group_by(var) %>%
  dplyr::summarise(mean_TP = mean(TP, na.rm = TRUE))

fit_model <- lm(mean_TP ~ poly(var, 2, raw = TRUE), data = df_mean)
data.pred <- data.frame( mean_TP = predict(fit_model, newdata = data.frame(var = df1$var)), var = df1$var)

df1$group <- ifelse(
  df1$var == 0.392307461244186| 
    df1$var >= 0.550336279069767,
  "sig", "no" )


fig1c <- ggplot() + 
  geom_point(data = df1, aes(x = var, y = TP, color = group ), size = 2, alpha = 0.08, show.legend = F ) +  
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "solid") +  
  labs(x = "Root diameter (mm)", y = "Standardized correlation effect") +
  theme(plot.background = element_blank(), panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() + 
  theme(axis.text = element_text(size = 11), axis.text.y = element_text(angle = 90, hjust = 0.5), axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12)) +
  theme(legend.position = c(0.53, 0.08), legend.title = element_blank(), legend.box.background = element_blank(), legend.background = element_blank(),
        legend.text = element_text(size = 12), legend.key = element_blank()) + 
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  geom_line(data = data.pred, aes(x = var, y = mean_TP, colour = "TP ~ Root length colonization"), size = 0.8, linetype = "solid",  show.legend = T) +  
  scale_color_manual(name = "", values = c("TP ~ Root length colonization" = "black","no" = "#386CB0", "sig" = "#ff8a00"), 
                     breaks = c("TP ~ Root length colonization"), labels = c("TP ~ Root length colonization" = "TP ~ Root length colonization")  )

fig1c


pp1 <- cowplot::plot_grid( fig1a,  fig1b,  fig1c,  
                           ncol = 3, nrow = 1, labels = c('(a)', "(b)", "(c)"), 
                           label_fontface = "bold", label_x = -0.03, label_y = 0.998 , align = c("hv"))
pp1

#1_fig1_12*4
fig1 <- cowplot::plot_grid( pp1, legend, align = "", axis = "tpl", rel_widths = c( 1.6, 0.28 )) 
fig1


figs2a <- ggplot(df1, aes( x=as.factor(round(var, 3)),y= TP )) + 
  geom_boxplot(width=0.5,outlier.shape = 1, outlier.size = .3,)+ geom_hline(yintercept=0)+
  scale_color_manual(name = "", values = c("no" = "#386CB0", "sig" = "#ff8a00"))+
  geom_jitter(aes(color = group ), shape=1, size=.3,position=position_jitter(0.05),alpha=0.2, show.legend = F)+ 
  xlab("Root diameter (mm)")+ylab("")+ ggtitle("Standardized bootstrapped coefficients of correlation between TP and root length colonization") + theme_bw()+theme(panel.grid=element_blank())+theme(plot.background = element_blank())+
  theme(axis.text.x = element_text(size = 9,color = "black",angle = 60, hjust = 1),axis.ticks.x = element_blank(),title = element_text(size = 11)) +
  theme(axis.text.y = element_text(size = 9,color = "black"))+theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 11)) +
  scale_x_discrete(labels = function(x) format(as.numeric(x), nsmall = 3))+
  annotate("text", 
           x = as.factor(round(unique(c(
             df1$var[df1$var == 0.392307461244186],
             df1$var[df1$var >= 0.550336279069767]
           )), 3)), 
           y = -0.35,  # 星号的位置在y轴下方
           label = "*", 
           size = 6, 
           color = "black")

figs2a

######## Distribution and sharing patterns of AM fungal OTUs ---------------
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
data = read.csv("field_data_otu_ra.csv", header=TRUE, row.names=1) 

otu_counts <- colSums(data > 0) 

df_plant <- data.frame( 
  Plant = names(otu_counts), 
  OTU_Count = as.numeric(otu_counts)
) 

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
  Num_OTUs = as.numeric(otu_dist_table)
)


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


#2_figs4_8*14
figs4 <- cowplot::plot_grid( gg_plant_otu,  gg_otu_dist , 
                             ncol = 1, nrow = 2, labels = c('(a)', "(b)", "(c)",'(d)', "(e)", "(f)", "(g)", "(h)", 
                                                            "(i)", "(j)","(k)", "(l)", "(m)","(n)","(o)", "(p)", 
                                                            "(q)","(r)","(s)", "(t)"), 
                             label_fontface = "bold", label_x = 0 , label_y = 1.01 , align = c("hv"))
figs4




######## Relationship between plant and AM fungal community composition ---------------
library(vegan)
library(picante)
library(ggtree)
library(ape)
library(phytools)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys")
data = read.csv("field_data_all.csv",header=TRUE, row.names=1) 
otu_rrarefy <- read.csv("field_data_otu_rrarefy.csv",header=TRUE, row.names=1)

otu_h <- decostand(otu_rrarefy, method = 'hellinger')    
vare.dist <- vegdist(otu_h, method="bray")

adonis2(vare.dist ~ data$Plant.species, permutations = 999) 

# Relationship between plant phylogeny and AM fungal communities
tree_raw <- read.tree("2021mqamf_35sp.nwk")
plot(tree_raw)

rooted_tree1 <- root(tree_raw, outgroup = "Juncus", node = NULL, resolve.root = T, interactive = FALSE, edgelabel = FALSE)
phylo  <- drop.tip(rooted_tree1, "Juncus")
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


######## Relationships between AM fungal OTUs and root diameter ---------------
library(ggplot2)
library(dplyr)
library(ggrepel)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
otu = read.csv("field_data_otu_rrarefy.csv", header=TRUE, row.names=1) 
data = read.csv("field_data_all.csv", header=TRUE, row.names=1) 

otu_rel <- sweep(otu, 1, rowSums(otu), FUN = "/") * 100
mean_rel <- colMeans(otu_rel)
otu_sel <- otu_rel[, mean_rel > 1]
total_RA_per_sample <- rowSums(otu_sel)
mean(total_RA_per_sample)

root_dia <- as.numeric(data$Root.diameter)

corr_values <- apply(otu_sel, 2, function(x) cor(x, root_dia, method = "pearson"))
p_values <- apply(otu_sel, 2, function(x) cor.test(x, root_dia, method = "pearson")$p.value)

corr_df <- data.frame(
  OTU = names(corr_values),
  r = corr_values,
  p = p_values,
  Mean_RA = mean_rel[names(corr_values)]
)

corr_df$significant <- ifelse(corr_df$p < 0.05, "Significant", "Not Significant")

corr_df$color <- ifelse(corr_df$r < 0 & corr_df$significant == "Significant", "Decrease",
                        ifelse(corr_df$r > 0 & corr_df$significant == "Significant", "Increase",
                               "Unchange"))

group =  read.csv("field_data_otu_taxon.csv", header=TRUE, row.names=1)
corr_df <- merge(corr_df, group, by.x = "OTU", by.y = "row.names", all.x = TRUE)
corr_df$label <- ifelse(corr_df$significant == "Significant" & corr_df$Mean_RA > 1, corr_df$OTU.id, NA)

#2_figs5_6*7
figs5 <- ggplot(corr_df, aes(x = r, y = -log10(p), size = Mean_RA, color = color)) +
  geom_point(alpha = 1, shape = 16) +
  scale_color_manual(values = c("Decrease" = "steelblue", "Increase" = "tomato", "Unchange" = "gray77")) +
  geom_text_repel(aes(label = label, size = Mean_RA/1.5, color = color),
                  max.overlaps = 5, fontface = "plain", na.rm = TRUE, show.legend = FALSE) +
  labs(x = "Pearson correlation coefficient", 
       y = "-Log10 (P value)", 
       size = "RA (%)",
       color = "OTUs") +
  scale_size_continuous(breaks = c(2, 4, 6)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.title = element_text(hjust = 0, size = 12),
        legend.text = element_text(size = 11),legend.key = element_blank() ,
        legend.position = c(0.12, 0.2),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(size = guide_legend(order = 1),
         color = guide_legend(order = 2))

figs5


######## AM fungal co-occurrence network advantage modules ---------------
library(circlize) 
library(dplyr)
library(readr)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
otu_anno <- read.csv("field_data_otu_network_module.csv",header=TRUE, row.names=1) 
edges <- read.csv("1_network_edge.csv",header=TRUE, row.names=1) 

otu_anno <- otu_anno %>%
  mutate(Module = ifelse(is.na(Module), "Other", as.character(Module)))  

otu_ids <- unique(c(edges$From, edges$To))
otu_anno <- otu_anno %>% filter(OTU_ID %in% otu_ids)
edges <- edges %>% filter(From %in% otu_ids & To %in% otu_ids)

module_colors <- c(
  "Module#1" = "#A83248",
  "Module#2" = "#BFA87F", 
  "Module#3" = "#55948D", 
  "Module#4" = "#377EB8", 
  "Module#5" = "#984EA3", 
  "Other" = "grey66"  
)


otu_anno <- otu_anno %>% 
  mutate(Module_label = ifelse(Module %in% c("Module#1", "Module#2", "Module#3", "Module#4", "Module#5"), Module, "Other"))

grid.col <- module_colors[otu_anno$Module_label]
names(grid.col) <- otu_anno$OTU_ID

#2_figs6_8*8
chordDiagram(
  x = edges,
  grid.col = grid.col,
  transparency = 0.1,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.4)
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector_name <- get.cell.meta.data("sector.index")
    label <- otu_anno$Annotation[otu_anno$OTU_ID == sector_name]
    label_col <- grid.col[sector_name]
    
    circos.text(
      x = CELL_META$xcenter,
      y = CELL_META$yplot[1] - mm_y(22), 
      labels = label,
      facing = "clockwise",
      niceFacing = TRUE,
      cex = 0.7,
      adj = c(0, 0.5),
      col = label_col  
    )
  },
  bg.border = NA
)

legend(
  x = "top",                 
  legend = names(module_colors),
  fill = module_colors,
  border = NA,
  bty = "n",
  cex = 0.9,
  title = NA,
  horiz = TRUE,              
  x.intersp = 1.1
)

######## Relationship between plant traits and AM fungal Richness & Module #1 ---------------
library(ggplot2)
library(gridExtra)
library(reshape)
library(gam)
library(segmented)
library(chngpt)
library(ggpmisc)
library(car)
library(doBy)
library(dplyr) 
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
data = read.csv("field_data_all.csv",header=TRUE, row.names=1) 

cols <- c("Amaryllidaceae"= "#6C8FA6", "Cyperaceae" = "#5F97C6", "Poaceae" = "#8FB9D1", 
          "Ranunculaceae" = "#BFA87F","Onagraceae" = "#B86A2E", "Brassicaceae"= "#B09AA8",  
          "Rosaceae" =  "#E0A16F", "Fabaceae" = "#FFA64D" , "Caryophyllaceae" = "#7FC1B5" , 
          "Polygonaceae"= "#6BA89F",  "Gentianaceae"= "#9e8c93","Lamiaceae"  = "#D94F30",
          "Plantaginaceae" = "#A83248", "Apiaceae" = "#B85A78", "Asteraceae" = "#E9746A") 

se <- function(x) sd(x)/sqrt(length(x))
conf_95 <- function(x) t.test(x, conf.level = 0.95)$conf.int[2]

outlier_test = outlierTest(lm( OTU.Richness ~ Root.diameter , data = data ))                         
outlier_test 
sumdata <- summaryBy( Root.diameter + OTU.Richness ~ Plant.species, data, FUN = c(mean, sd, se))
sumdata
sumdata$Family <- data$Plant.family[match(sumdata$Plant.species, data$Plant.species)]
sumdata$Root.diameter <- sumdata$Root.diameter.mean
sumdata$OTU.Richness <- sumdata$OTU.Richness.mean
colnames(sumdata)

fig2a <- ggplot(sumdata, aes(y = OTU.Richness , x = Root.diameter,  color = Family  ))+      
  scale_color_manual(values = cols)+  labs(color="Plant family") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "OTU richness" )+
  geom_errorbarh(aes(xmin = Root.diameter + Root.diameter.se, xmax = Root.diameter - Root.diameter.se), height=0,  show.legend = F) +
  geom_errorbar( aes(ymin = OTU.Richness + OTU.Richness.se, ymax = OTU.Richness - OTU.Richness.se), width =0,  show.legend = F) +
  geom_point(size = 2.5, alpha = 1, shape = 16, aes(color = Family), show.legend = F ) + 
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  stat_poly_eq(  aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")), data = data, formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.98,label.y = 0.02, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  theme_bw() + theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5),
                     legend.text=element_text(size=11), legend.title = element_text(size =12, hjust=0.1),  axis.title=element_text(size=12)) + 
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", label.hjust = 1))

fig2a


outlier_test = outlierTest(lm( P.concentrations ~ OTU.Richness , data = data)) 
outlier_test 
sumdata <- summaryBy( P.concentrations + OTU.Richness ~ Plant.species, data, FUN = c(mean, sd, se ))
sumdata
sumdata$Family <- data$Plant.family[match(sumdata$Plant.species, data$Plant.species)]
sumdata$P.concentrations <- sumdata$P.concentrations.mean
sumdata$OTU.Richness <- sumdata$OTU.Richness.mean
colnames(sumdata)

fig2c <- ggplot(sumdata, aes(y = P.concentrations, x = OTU.Richness, color = Family  ))+      
  #geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "dashed" , data = data )+ 
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")), data = data, formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.98,label.y = 0.98, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  scale_color_manual(values = cols)+  labs(color="Plant family")+
  scale_x_continuous(name ="OTU richness") + scale_y_continuous(name = "TP (mg/g)" )+
  geom_errorbarh(aes(xmin = OTU.Richness + OTU.Richness.se, xmax = OTU.Richness - OTU.Richness.se),  height=0, show.legend = F) +
  geom_errorbar(aes(ymin = P.concentrations + P.concentrations.se, ymax = P.concentrations - P.concentrations.se), width=0, show.legend = F) +
  geom_point(size =  2.5, alpha = 1, shape = 16, aes(color = Family), show.legend = F ) + 
  theme_bw() + theme(axis.text=element_text(size=11),  axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=11), legend.title = element_text(size=12, hjust=0.1),  axis.title=element_text(size=12)) + #
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", label.hjust = 1))
fig2c


outlier_test = outlierTest(lm( Family.Richness ~ Root.diameter , data = data)) 
outlier_test 
sumdata <- summaryBy( Root.diameter + Family.Richness ~ Plant.species, data, FUN = c(mean, sd, se))
sumdata
sumdata$Family <- data$Plant.family[match(sumdata$Plant.species, data$Plant.species)] 
sumdata$Root.diameter <- sumdata$Root.diameter.mean
sumdata$Family.Richness <- sumdata$Family.Richness.mean
colnames(sumdata)

fig2b <- ggplot(sumdata, aes(y = Family.Richness , x = Root.diameter,  color = Family  ))+      
  #geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "dashed" , data = data)+ 
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")),  data = data, formula = y ~ poly(x, 1, raw = TRUE), parse = TRUE,label.x = 0.98,label.y = 0.98, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  scale_color_manual(values = cols)+  labs(color="Plant family") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "Family richness" )+
  geom_errorbarh(aes(xmin = Root.diameter + Root.diameter.se, xmax = Root.diameter - Root.diameter.se),  height = 0, show.legend = F) +
  geom_errorbar(aes(ymin = Family.Richness + Family.Richness.se, ymax = Family.Richness - Family.Richness.se), width=0, show.legend = F) +
  geom_point(size =  2.5, alpha = 1, shape = 16, aes(color = Family), show.legend = F ) + 
  theme_bw() + theme(axis.text=element_text(size = 11), axis.text.y = element_text(angle = 90, hjust = 0.5),
                     legend.text=element_text(size= 11 ), legend.title = element_text(hjust=0.1),  axis.title=element_text(size=12)) + #
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", label.hjust = 0))
fig2b


outlier_test = outlierTest(lm( P.concentrations ~ Family.Richness , data = data )) 
outlier_test 
sumdata <- summaryBy( P.concentrations + Family.Richness ~ Plant.species, data, FUN = c(mean, sd, se))
sumdata
sumdata$Family <- data$Plant.family[match(sumdata$Plant.species, data$Plant.species)] 
sumdata$P.concentrations <- sumdata$P.concentrations.mean
sumdata$Family.Richness <- sumdata$Family.Richness.mean
colnames(sumdata)

fig2d <- ggplot(sumdata, aes(y = P.concentrations , x = Family.Richness, color = Family  )) +      
  scale_color_manual(values = cols )+  labs(color="Plant family")+
  scale_x_continuous(name ="Family richness") + scale_y_continuous(name = "TP (mg/g)" ) +
  geom_errorbarh(aes(xmin = Family.Richness + Family.Richness.se, xmax = Family.Richness - Family.Richness.se), height=0, show.legend = F) +
  geom_errorbar(aes(ymin = P.concentrations + P.concentrations.se, ymax = P.concentrations - P.concentrations.se), width=0, show.legend = F ) +
  geom_point(size =  2.5, alpha = 1, shape = 16, aes(color = Family), show.legend = F ) + 
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ # #3366FF #C53F1C
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")),  data = data, formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.98,label.y = 0.98, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  theme_bw() + theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=11), legend.title = element_text(size=12, hjust=0.15),  axis.title=element_text(size=12)) + #
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", label.hjust = 0))
fig2d


data_s = data[-c(38, 54, 35), ]
outlier_test = outlierTest(lm( RA.of.Module..1 ~ Root.diameter , data = data_s )) 
outlier_test 
sumdata <- summaryBy( Root.diameter + RA.of.Module..1 ~ Plant.species, data_s , FUN = c(mean, sd, se))
sumdata
sumdata$Family <- data$Plant.family[match(sumdata$Plant.species, data$Plant.species)] 
sumdata$Root.diameter <- sumdata$Root.diameter.mean
sumdata$RA.of.Module..1 <- sumdata$RA.of.Module..1.mean
colnames(sumdata)

fig2e <- ggplot(sumdata, aes(y = RA.of.Module..1 , x = Root.diameter,  color = Family  )) +      
  scale_color_manual(values = cols)+  labs(color="Plant family") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "RA of Module #1 (%)" )+
  geom_errorbarh(aes(xmin = Root.diameter + Root.diameter.se, xmax = Root.diameter - Root.diameter.se), height =0,  show.legend = F) +
  geom_errorbar(aes(ymin = RA.of.Module..1 + RA.of.Module..1.se, ymax = RA.of.Module..1 - RA.of.Module..1.se), width = 0,  show.legend = F) +
  geom_point(size =  2.5, alpha = 1, shape = 16, aes(color = Family), show.legend = F ) + 
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data_s )+ # #3366FF #C53F1C
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")), data = data_s , formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.98,label.y = 0.02, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  theme_bw() + theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), legend.text=element_text(size=11), legend.title = element_text(size=12, hjust=0.8),  axis.title=element_text(size=12)) + #
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))
fig2e


outlier_test = outlierTest(lm( P.concentrations ~ RA.of.Module..1 , data = data_s )) 
outlier_test 
sumdata <- summaryBy( P.concentrations + RA.of.Module..1 ~ Plant.species, data_s, FUN = c(mean, sd, se))
sumdata
sumdata$Family <- data$Plant.family[match(sumdata$Plant.species, data$Plant.species)] 
sumdata$P.concentrations <- sumdata$P.concentrations.mean
sumdata$RA.of.Module..1 <- sumdata$RA.of.Module..1.mean
colnames(sumdata)

fig2f <- ggplot(sumdata, aes(y = P.concentrations , x = RA.of.Module..1, color = Family  ))+      
  scale_color_manual(values = cols, breaks = c("Asteraceae",  "Plantaginaceae", "Lamiaceae",  
                                               "Gentianaceae",   "Rosaceae", "Caryophyllaceae",
                                               "Ranunculaceae", "Poaceae","Cyperaceae", "Amaryllidaceae" ) ) +  labs(color="Plant family")+ 
  scale_x_continuous(name ="RA of Module #1 (%)") + scale_y_continuous(name = "TP (mg/g)" )+
  geom_errorbarh(aes(xmin = RA.of.Module..1 + RA.of.Module..1.se, xmax = RA.of.Module..1 - RA.of.Module..1.se),  show.legend = F) +
  geom_errorbar(aes(ymin = P.concentrations + P.concentrations.se, ymax = P.concentrations - P.concentrations.se),  show.legend = F ) +
  geom_point(size =  2.5, alpha = 1, shape = 16, aes(color = Family), show.legend = F ) + 
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data_s )+ # #3366FF #C53F1C
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")), data = data_s, formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.98,label.y = 0.02, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  theme_bw() + theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), legend.text=element_text(size=11), legend.title = element_text(size=12, hjust=0.15),  axis.title=element_text(size=12)) + #
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", label.hjust = 0)) #  + coord_cartesian(xlim = c(0, 45))
fig2f


legend <- cowplot::get_legend(fig2f) 

###
colnames(data)
dat <- data[ , c( 3, 5 , 6, 8)]
colnames(dat)
head(dat)

set.seed(12345)
bootari <- function(DF, model, var,
                    window = 85,
                    NsimB = 500,
                    ncore = max(1, parallel::detectCores() - 2)) {
  stopifnot(var %in% colnames(DF))
  xvar <- DF[[var]]
  if (is.data.frame(xvar)) xvar <- xvar[[1]]
  if (!is.numeric(xvar)) stop("var must be numeric")
  DF[[var]] <- xvar
  DF <- DF[order(DF[[var]]), ]
  vectorVAR <- sort(unique(DF[[var]]))
  nVAR <- length(vectorVAR)
  
  if (nVAR <= window) {
    stop("Not enough unique values in var for the chosen window size")
  }
  
  library(lme4)
  library(boot)
  library(doSNOW)
  
  cl <- makeCluster(ncore)
  registerDoSNOW(cl)
  
  bt <- NULL
  varNote <- numeric()
  
  for (i in 1:(nVAR - window)) {
    
    try({
      
      winVAR <- vectorVAR[i:(i + window)]
      
      dfi <- DF[DF[[var]] %in% winVAR, ]
      
      
      dfi$P.concentrations <-
        scale(dfi$P.concentrations, center = TRUE, scale = TRUE)
      dfi$Family.Richness <-
        scale(dfi$Family.Richness, center = TRUE, scale = TRUE)
      dfi$Root.diameter <- dfi$Root.diameter 
      
      fit <- lmer(as.formula(model), data = dfi, REML = FALSE)
      
      bt_i <- bootMer(
        fit,
        FUN = fixef,
        nsim = NsimB,
        parallel = "snow",
        cl = cl
      )$t
      
      
      bt <- rbind(bt, bt_i)
      varNote <- c(varNote, rep(mean(winVAR), NsimB))
      
    }, silent = TRUE)
  }
  
  stopCluster(cl)
  
  
  bt <- as.data.frame(bt)
  colnames(bt) <- names(fixef(fit))
  bt$var <- varNote
  
  return(bt)
}


model_1 <- "P.concentrations ~ Family.Richness + (1|Block)"
var <- "Root.diameter"

btaridity <- bootari(
  DF = dat,
  model = model_1,
  var = var,
  window = 85,
  NsimB = 500 )
colnames(btaridity) <- c("Intcp","TP","var")

model_11 <-"P.concentrations ~ Family.Richness + (1|Block)"
var1<-"Root.diameter"
btaridity1 <- bootari(
  DF = dat,
  model = model_11,
  var = var1,
  window = 85,
  NsimB = 500 )
colnames(btaridity1) <- c("Intcp","TP","var")

lsbt<-list(bt1=btaridity, bt2=btaridity1)
lsgood<-lsbt

for (n in 1:length(lsbt)){
  btn<-lsbt[[n]]
  btn<-btn[complete.cases(btn),]
  for (i in 1:length(unique(btn$var))){
    bti<-btn[btn$var==unique(btn$var)[i],]
    for(j in 1:(length(colnames(bti))-1)){
      x<-bti[,j]
      upci<-quantile(x, 0.95)
      lowci<-quantile(x,0.05)
      x[x>=upci | x<=lowci]<-NA
      bti[,j]<-x
    }
    btn[btn$var==unique(btn$var)[i],]<-bti
  }
  lsgood[[n]]<-btn
}

df1=lsbt[[1]]
colnames(df1)
length(unique(df1$var))
table(df1$var)

#write.csv(df1, "1_150_MW_85_TP_FR_RAD.csv", quote=F)

df1=read.csv("1_150_MW_85_TP_FR_RAD.csv",head=T,row.name=1)


df1$group <- ifelse(
  df1$var == 0.497104844965116 | 
    df1$var == 0.504669031011628 | 
    
    df1$var == 0.527472906976744 | 
    df1$var == 0.535004302325581 | 
    df1$var == 0.542504011627907 |  
    df1$var == 0.550336279069767 | 
    
    df1$var == 0.576177945732558 |   
    df1$var == 0.585205852709302 |   
    df1$var == 0.59427523255814 |   
    df1$var == 0.603415465116279 | 
    
    df1$var == 0.612975930232558 |   
    df1$var == 0.622808488372093 |   
    df1$var == 0.632571860465116 |   
    df1$var == 0.642643953488372 |   
    df1$var == 0.652953837209302 |  
    df1$var == 0.663310813953488 |  
    
    df1$var == 0.686733197674419 | 
    df1$var == 0.726169050383721 ,

  "no", "sig" )

figs2b <- ggplot(df1, aes( x=as.factor(round(var, 3)),y= TP )) + 
  geom_boxplot(width=0.5, outlier.shape = 1, outlier.size = .3,)+ geom_hline(yintercept=0)+
  scale_color_manual(name = "", values = c("no" = "#386CB0", "sig" = "#ff8a00"))+
  geom_jitter(aes(color = group ), shape=1, size=.3,position=position_jitter(0.05),alpha=0.2, show.legend = F)+ 
  xlab("Root diameter (mm)")+ylab("")+ ggtitle("Standardized bootstrapped coefficients of correlation between TP and family richness") + theme_bw()+theme(panel.grid=element_blank())+theme(plot.background = element_blank())+
  theme(axis.text.x = element_text(size = 9,color = "black",angle = 60, hjust = 1), axis.ticks.x = element_blank(),title = element_text(size = 11)) +
  theme(axis.text.y = element_text(size = 9,color = "black")) + theme(axis.title.x = element_text(size = 12)) +
  theme( axis.title.y = element_text(size = 11) ) +
  scale_x_discrete(labels = function(x) format(as.numeric(x), nsmall = 3))+
  annotate("text", 
           x = as.factor(round(unique(df1$var[!(df1$var %in% c(
             0.497104844965116, 0.504669031011628, 0.527472906976744, 0.535004302325581, 0.542504011627907, 0.550336279069767, 0.576177945732558,
             0.585205852709302, 0.59427523255814, 0.603415465116279, 0.612975930232558, 0.622808488372093, 0.632571860465116, 0.642643953488372,
             0.652953837209302, 0.663310813953488, 0.686733197674419, 0.726169050383721
           ))]), 3)), 
           y = -0.28,  
           label = "*", 
           size = 5, 
           color = "black")

figs2b


###
colnames(data)
data_s = data[-c(38, 54, 35), ]
dat <- data_s[ , c( 3, 5 , 6, 10)]
colnames(dat)
head(dat)

set.seed(12345)
bootari <- function(DF, model, var,
                    window = 85,
                    NsimB = 500,
                    ncore = max(1, parallel::detectCores() - 2)) {
  stopifnot(var %in% colnames(DF))
  xvar <- DF[[var]]
  if (is.data.frame(xvar)) xvar <- xvar[[1]]
  if (!is.numeric(xvar)) stop("var must be numeric")
  DF[[var]] <- xvar
  DF <- DF[order(DF[[var]]), ]
  vectorVAR <- sort(unique(DF[[var]]))
  nVAR <- length(vectorVAR)
  
  if (nVAR <= window) {
    stop("Not enough unique values in var for the chosen window size")
  }
  
  library(lme4)
  library(boot)
  library(doSNOW)
  
  cl <- makeCluster(ncore)
  registerDoSNOW(cl)
  
  bt <- NULL
  varNote <- numeric()
  
  for (i in 1:(nVAR - window)) {
    
    try({
      
      winVAR <- vectorVAR[i:(i + window)]
      
      dfi <- DF[DF[[var]] %in% winVAR, ]
      
      
      dfi$P.concentrations <-
        scale(dfi$P.concentrations, center = TRUE, scale = TRUE)
      dfi$RA.of.Module..1 <-
        scale(dfi$RA.of.Module..1, center = TRUE, scale = TRUE)
      dfi$Root.diameter <- dfi$Root.diameter 
      
      fit <- lmer(as.formula(model), data = dfi, REML = FALSE)
      
      bt_i <- bootMer(
        fit,
        FUN = fixef,
        nsim = NsimB,
        parallel = "snow",
        cl = cl
      )$t
      
      
      bt <- rbind(bt, bt_i)
      varNote <- c(varNote, rep(mean(winVAR), NsimB))
      
    }, silent = TRUE)
  }
  
  stopCluster(cl)
  
  
  bt <- as.data.frame(bt)
  colnames(bt) <- names(fixef(fit))
  bt$var <- varNote
  
  return(bt)
}


model_1 <- "P.concentrations ~ RA.of.Module..1 + (1|Block)"
var <- "Root.diameter"

btaridity <- bootari(
  DF = dat,
  model = model_1,
  var = var,
  window = 85,
  NsimB = 500 )
colnames(btaridity) <- c("Intcp","TP","var")

model_11 <-"P.concentrations ~ RA.of.Module..1 + (1|Block)"
var1<-"Root.diameter"
btaridity1 <- bootari(
  DF = dat,
  model = model_11,
  var = var1,
  window = 85,
  NsimB = 500 )
colnames(btaridity1) <- c("Intcp","TP","var")

lsbt<-list(bt1=btaridity, bt2=btaridity1)
lsgood<-lsbt

for (n in 1:length(lsbt)){
  btn<-lsbt[[n]]
  btn<-btn[complete.cases(btn),]
  for (i in 1:length(unique(btn$var))){
    bti<-btn[btn$var==unique(btn$var)[i],]
    for(j in 1:(length(colnames(bti))-1)){
      x<-bti[,j]
      upci<-quantile(x, 0.95)
      lowci<-quantile(x,0.05)
      x[x>=upci | x<=lowci]<-NA
      bti[,j]<-x
    }
    btn[btn$var==unique(btn$var)[i],]<-bti
  }
  lsgood[[n]]<-btn
}

df1=lsbt[[1]]
colnames(df1)
length(unique(df1$var))
table(df1$var)

#write.csv(df1, "1_147_MW_85_TP_M1_RAD.csv", quote=F)

df1=read.csv("1_147_MW_85_TP_M1_RAD.csv",head=T,row.name=1) 

df_mean <- df1 %>% 
  dplyr::group_by(var) %>% 
  dplyr::summarise(mean_TP = mean(TP, na.rm = TRUE)) 


fit_model <- lm(mean_TP ~ poly(var, 2, raw = TRUE), data = df_mean)
data.pred <- data.frame( mean_TP = predict(fit_model, newdata = data.frame(var = df1$var)), var = df1$var)

df1$group <- ifelse(
  df1$var >= 0.607990755813953,
  "sig", "no" )

fig2g <- ggplot() + 
  geom_point(data = df1, aes(x = var, y = TP, color = group ), size = 2.5, alpha = 0.08, show.legend = F ) +  
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "solid") +  
  labs(x = "Root diameter (mm)", y = "Standardized correlation effect") +
  theme(plot.background = element_blank(), panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() + 
  theme(axis.text = element_text(size = 11), axis.text.y = element_text(angle = 90, hjust = 0.5), axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12)) +
  theme(legend.position = c(0.36,0.95), legend.title = element_blank(), legend.box.background = element_blank(), legend.background = element_blank(),
        legend.text = element_text(size = 12), legend.key = element_blank()) +
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  geom_line(data = data.pred, aes(x = var, y = mean_TP, colour = "TP ~ RA of Module #1"), size = 0.8, linetype = "solid",  show.legend = T) +  
  scale_color_manual(name = "", values = c("TP ~ RA of Module #1" = "black","no" = "#386CB0", "sig" = "#ff8a00"),
                     breaks = c("TP ~ RA of Module #1"),
                     labels = c("TP ~ RA of Module #1" = "TP ~ RA of Module #1")  )

fig2g


figs2c <- ggplot(df1, aes( x=as.factor(round(var, 3)),y= TP )) + 
  geom_boxplot(width=0.5,outlier.shape = 1,outlier.size = .3,)+ geom_hline(yintercept=0)+
  scale_color_manual(name = "", values = c("no" = "#386CB0", "sig" = "#ff8a00"))+
  geom_jitter(aes(color = group ), shape=1, size=.3, position=position_jitter(0.05), alpha=0.2, show.legend = F)+ 
  xlab("Root diameter (mm)")+ylab("")+ ggtitle("Standardized bootstrapped coefficients of correlation between TP and RA of Module #1") + theme_bw()+theme(panel.grid=element_blank())+theme(plot.background = element_blank())+
  theme(axis.text.x = element_text(size = 9,color = "black",angle = 60, hjust = 1),axis.ticks.x = element_blank(),title = element_text(size = 11)) +
  theme(axis.text.y = element_text(size = 9,color = "black"))+theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 11)) +
  scale_x_discrete(labels = function(x) format(as.numeric(x), nsmall = 3)) +
  annotate("text", 
           x = as.factor(round(unique(c(
             df1$var[df1$var >= 0.607990755813953]
           )), 3)), 
           y = -0.35,  
           label = "*", 
           size = 6, 
           color = "black")
figs2c 

#2_fig.s2_12*8
figs2 <- cowplot::plot_grid( figs2a,  figs2b, figs2c, 
                           ncol = 1, nrow = 3, labels = c( "(a)", "(b)","(c)",'(d)', "(e)", "(f)", "(g)", "(h)", 
                                                           "(i)", "(j)","(k)", "(l)", "(m)", "(n)", "(o)", "(p)",  
                                                           "(q)","(r)","(s)", "(t)"), 
                           label_fontface = "bold", label_x = 0.01, label_y = 1.0 , align = c("hv"))
figs2


fig2.1 <- cowplot::plot_grid( fig2a,  fig2b,  fig2c, fig2d, 
                           ncol = 4, nrow = 1, labels = c(  "(a)", "(b)","(c)",'(d)', "(e)", "(f)", "(g)" ),
                           label_fontface = "bold", label_x = -0.02, label_y = 0.998 , align = c("hv"))
fig2.1

fig2.2 <- cowplot::plot_grid( fig2e,  fig2f,  fig2g, legend,
                              ncol = 4, nrow = 1, labels = c(  "(e)", "(f)", "(g)" ),
                              label_fontface = "bold", label_x = -0.02, label_y = 0.998 , align = c("hv"))
fig2.2

#1_fig2_14*8
fig2 <- cowplot::plot_grid(
  fig2.1, 
  fig2.2, 
  ncol = 1,          
  align = "hv",       
  axis = "tb",      
  rel_heights = c(1, 1) )

fig2

######## Phylogenetic generalized least squares regressions ----------------
library(car)
library(caper) 
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
data = read.csv("field_data_mean.csv", header=TRUE, row.names=1) 

tree_raw <- read.tree("2021mqamf_35sp.nwk") 
plot(tree_raw)
rooted_tree1 <- root(tree_raw, outgroup = "Juncus", node = NULL, resolve.root = T, interactive = FALSE, edgelabel = FALSE)
plot(rooted_tree1)
rooted_tree  <- drop.tip(rooted_tree1, "Juncus")
plot(rooted_tree)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data , names.col = "Plant.species", vcv = TRUE, vcv.dim = 2)

colnames(data)
###
outlier_test = outlierTest( lm( Root.length.colonization ~ Root.diameter , data =  data)) 
outlier_test
pgls_model <- pgls(  Root.length.colonization ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)

outlier_test = outlierTest(lm( P.concentrations ~ Root.length.colonization, data =  data)) 
outlier_test
pgls_model <- pgls(P.concentrations ~ Root.length.colonization, data = comp_data, lambda = "ML") 
summary(pgls_model)

###
outlier_test = outlierTest(lm( OTU.Richness ~ Root.diameter , data =  data)) 
outlier_test
pgls_model <- pgls(  OTU.Richness ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)

outlier_test = outlierTest(lm( P.concentrations ~ OTU.Richness, data =  data)) 
outlier_test
pgls_model <- pgls( P.concentrations ~ OTU.Richness , data = comp_data, lambda = "ML") 
summary(pgls_model)

###
outlier_test = outlierTest(lm( Family.Richness ~ Root.diameter , data =  data)) 
outlier_test
pgls_model <- pgls( Family.Richness ~ Root.diameter , data = comp_data, lambda = "ML")
summary(pgls_model)

outlier_test = outlierTest(lm( P.concentrations ~ Family.Richness, data =  data)) 
outlier_test
pgls_model <- pgls( P.concentrations ~ Family.Richness , data = comp_data, lambda = "ML")
summary(pgls_model)

###
outlier_test = outlierTest(lm( RA.of.Module..1 ~ Root.diameter , data =  data)) 
outlier_test
pgls_model <- pgls( RA.of.Module..1 ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)

outlier_test = outlierTest(lm( P.concentrations ~ RA.of.Module..1, data =  data)) 
outlier_test
pgls_model <- pgls( P.concentrations ~ RA.of.Module..1 , data = comp_data, lambda = "ML") 
summary(pgls_model)
