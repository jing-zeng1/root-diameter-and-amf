

######## Effects of plant species identity on plant traits and AM fungal variables ----------------
library(nlme)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/greenhouse experiments")  
data1 = read.csv("greenhouse_data_1.csv",header=TRUE, row.names=1)
data2 = read.csv("greenhouse_data_2.csv",header=TRUE, row.names=1)

###
fm <- lme(Root.diameter ~ Plant.species * Inoculation.treatment, data = data1 , random=~1|Block)
anova(fm)
fm <- lme(Root.diameter ~ Plant.species , data = data2 , random=~1|Block) 
anova(fm)
fm <- lme(Root.diameter_NM ~ Plant.species , data = data2 , random=~1|Block)
anova(fm)

###
fm <- lme(P.concentrations ~ Plant.species * Inoculation.treatment , data = data1 , random=~1|Block)
anova(fm)
fm <- lme( P.concentrations ~ Plant.species , data = data2 , random=~1|Block)
anova(fm)
fm <- lme( P.concentrations_NM ~ Plant.species , data = data2 , random=~1|Block)
anova(fm)

###
fm <- lme(Total.biomass ~ Plant.species * Inoculation.treatment , data = data1 , random=~1|Block) 
anova(fm)
fm <- lme( Total.biomass ~ Plant.species , data = data2 , random=~1|Block)
anova(fm)
fm <- lme( Total.biomass_NM ~ Plant.species , data = data2 , random=~1|Block)
anova(fm)

###
fm <- lme( Root.length.colonization ~ Plant.species  , data = data2 , random=~1|Block )
anova(fm)

fm <- lme( MPR ~ Plant.species  , data = data2 , random=~1|Block)
anova(fm)
fm <- lme( MGR ~ Plant.species  , data = data2 , random=~1|Block)
anova(fm)
fm <- lme( MRR ~ Plant.species  , data = data2 , random=~1|Block)
anova(fm)


######## Plant phylogenetic signals of plant traits and AM fungal variables ----------------
library(plyr) 
library(ggtree)
library(ape)
library(picante)
library(phylosignal)
library(phylobase)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/greenhouse experiments")  
data = read.csv("greenhouse_data_2.csv",header=TRUE, row.names=1)
colnames(data)

data_mean <- ddply(data, c("Plant.species"),summarise,   
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
tree_raw <- read.tree("2022myco_45sp.nwk") 
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
#write.csv(phy.sig, '1_greenhouse_data_physig.csv')  

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
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/greenhouse experiments")  
data = read.csv("greenhouse_data_2.csv",header=TRUE, row.names=1)

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
sumdata$Plant.family <- data$Plant.family[match(sumdata$Plant.species, data$Plant.species)]
sumdata$Root.diameter <- sumdata$Root.diameter.mean
sumdata$Root.length.colonization <- sumdata$Root.length.colonization.mean
colnames(sumdata)

figs7a <- ggplot(sumdata, aes(y = Root.length.colonization , x = Root.diameter, color = Plant.family  ))+      
  scale_color_manual(values = cols ) +  labs(color="Plant family")+
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "Root length colonization (%)" )+
  geom_errorbarh(aes(xmin = Root.diameter + Root.diameter.se, xmax = Root.diameter - Root.diameter.se), height = 0, show.legend = F) + 
  geom_errorbar(aes(ymin = Root.length.colonization + Root.length.colonization.se, ymax = Root.length.colonization - Root.length.colonization.se), width = 0, show.legend = F) +
  geom_point(size =  2.5, alpha = 1, shape = 16, aes(color = Plant.family), show.legend = F ) + 
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")), data = data , formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.02,label.y = 0.98, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  theme_bw() + theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=11), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))
figs7a

data_s = data[-c(221,224), ]
outlier_test = outlierTest(lm( P.concentrations ~ Root.length.colonization , data = data_s)) 
outlier_test 
sumdata <- summaryBy( P.concentrations + Root.length.colonization ~ Plant.species, data_s, FUN = c(mean, sd, se ))
sumdata
sumdata$Plant.family <- data_s$Plant.family[match(sumdata$Plant.species, data_s$Plant.species)]
sumdata$P.concentrations <- sumdata$P.concentrations.mean
sumdata$Root.length.colonization <- sumdata$Root.length.colonization.mean
colnames(sumdata)

figs7b <- ggplot(sumdata, aes(y = P.concentrations , x = Root.length.colonization, color = Plant.family  ))+      
  scale_x_continuous(name ="Root length colonization (%)") + scale_y_continuous(name = "TP (mg/g)" )+
  geom_errorbarh(aes(xmin = Root.length.colonization + Root.length.colonization.se, xmax = Root.length.colonization - Root.length.colonization.se), height = 0,show.legend = F) +
  geom_errorbar(aes(ymin = P.concentrations + P.concentrations.se, ymax = P.concentrations - P.concentrations.se), width = 0, show.legend = F ) +
  geom_point(size =  2.5, alpha = 1, shape = 16, aes(color = Plant.family), show.legend = F ) + 
  scale_color_manual(values = cols, breaks = c("Asteraceae",  "Apiaceae", "Plantaginaceae", 
                                               "Caryophyllaceae", "Polygonaceae", "Fabaceae", "Onagraceae", "Rosaceae",
                                               "Brassicaceae", "Ranunculaceae", "Poaceae", "Cyperaceae"  ) )+  labs(color="Plant family")+
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data_s )+ 
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")), data = data_s, formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.02, label.y = 0.98, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  theme_bw() + theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=11), legend.title = element_text(size=12, hjust=0.15),  axis.title=element_text(size=12)) + #
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", label.hjust = 0))
figs7b 


data_s = data[-c(174,35,171,173,175,172), ]
outlier_test = outlierTest(lm( Total.biomass ~ Root.length.colonization , data = data_s)) 
outlier_test 
sumdata <- summaryBy( Total.biomass + Root.length.colonization ~ Plant.species, data_s, FUN = c(mean, sd, se ))
sumdata
sumdata$Plant.family <- data_s$Plant.family[match(sumdata$Plant.species, data_s$Plant.species)]
sumdata$Total.biomass <- sumdata$Total.biomass.mean
sumdata$Root.length.colonization <- sumdata$Root.length.colonization.mean
colnames(sumdata)

figs7c <- ggplot(sumdata, aes(y = Total.biomass , x = Root.length.colonization, color = Plant.family  ))+      
  scale_x_continuous(name ="Root length colonization (%)") + scale_y_continuous(name = "TB (g/plant)" )+
  geom_errorbarh(aes(xmin = Root.length.colonization + Root.length.colonization.se, xmax = Root.length.colonization - Root.length.colonization.se), height = 0,show.legend = F) +
  geom_errorbar(aes(ymin = Total.biomass + Total.biomass.se, ymax = Total.biomass - Total.biomass.se), width = 0, show.legend = F ) +
  geom_point(size =  2.5, alpha = 1, shape = 16, aes(color = Plant.family), show.legend = F ) + 
  scale_color_manual(values = cols, breaks = c("Asteraceae",  "Plantaginaceae", "Lamiaceae",
                                               "Gentianaceae",   "Rosaceae", "Caryophyllaceae",
                                               "Ranunculaceae", "Poaceae","Cyperaceae", "Amaryllidaceae" ) )+  labs(color="Plant family")+
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data_s )+ 
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")), data = data_s, formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.02, label.y = 0.98, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  theme_bw() + theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=11), legend.title = element_text(size=12, hjust=0.15),  axis.title=element_text(size=12)) + #
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", label.hjust = 0))
figs7c 

legend <- cowplot::get_legend(figs7b)

###
data_s = data[-c(221,224), ]
colnames(data_s)
dat <- data_s[ , c( 4 , 5, 7, 11)]
colnames(dat)
head(dat)

set.seed(12345)
bootari <- function(DF, model, var,
                    window = 115,
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
  window = 115,
  NsimB = 500 )
colnames(btaridity) <- c("Intcp","TP","var")

model_11 <-"P.concentrations ~ Root.length.colonization+(1|Block)"
var1<-"Root.diameter"
btaridity1 <- bootari(
  DF = dat,
  model = model_11,
  var = var1,
  window = 115,
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

#write.csv(df1, "1_223_MW_115_TP_RLC_RAD.csv", quote=F)

df1=read.csv("1_223_MW_115_TP_RLC_RAD.csv",head=T,row.name=1)

df_mean <- df1 %>%
  dplyr::group_by(var) %>%
  dplyr::summarise(mean_TP = mean(TP, na.rm = TRUE))

fit_model <- lm(mean_TP ~ poly(var, 2, raw = TRUE), data = df_mean)
data.pred <- data.frame( mean_TP = predict(fit_model, newdata = data.frame(var = df1$var)), var = df1$var)

df1$group <- ifelse(
  
    df1$var == 0.244455172413793 | 
    df1$var == 0.249219827586207 | 
    df1$var == 0.253882327586207 | 
    df1$var == 0.255735775862069 | 
    df1$var == 0.257581896551724 | 
    
    df1$var >= 0.259435344827586, 
  "sig" , "no" )


figs7d <- ggplot() + 
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

figs7d

figs7f <- ggplot(df1, aes( x=as.factor(round(var, 3)),y= TP )) + 
  geom_boxplot(width=0.5,outlier.shape = 1, outlier.size = .3,)+ geom_hline(yintercept=0)+
  scale_color_manual(name = "", values = c("no" = "#386CB0", "sig" = "#ff8a00"))+
  geom_jitter(aes(color = group ), shape=1, size=.3,position=position_jitter(0.05),alpha=0.2, show.legend = F)+ 
  xlab("Root diameter (mm)")+ylab("")+ ggtitle("Standardized bootstrapped coefficients of correlation between TP and root length colonization") + theme_bw()+theme(panel.grid=element_blank())+theme(plot.background = element_blank())+
  theme(axis.text.x = element_text(size = 9,color = "black",angle = 60, hjust = 1),axis.ticks.x = element_blank(),title = element_text(size = 11)) +
  theme(axis.text.y = element_text(size = 9,color = "black"))+theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 11)) +
  scale_x_discrete(labels = function(x) format(as.numeric(x), nsmall = 3))+
  annotate(
    "text", 
    x = as.factor(
      round(
        unique(
          df1$var[
            
              df1$var == 0.244455172413793 | 
                df1$var == 0.249219827586207 | 
                df1$var == 0.253882327586207 | 
                df1$var == 0.255735775862069 | 
                df1$var == 0.257581896551724 | 
                
                df1$var >= 0.259435344827586
              
          ]  ), 3 )),  y = -0.32,   label = "*", size = 5, color = "black" )

figs7f


###
data_s <- data[-c(174,35,171,173,175,172), ]
colnames(data)
dat <- data_s[ , c(  4 , 5, 9, 11 )]
colnames(dat)
head(dat)

set.seed(12345)
bootari <- function(DF, model, var,
                    window = 115,
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

      dfi$Total.biomass <-
        scale(dfi$Total.biomass, center = TRUE, scale = TRUE)
      
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


model_1 <- "Total.biomass ~ Root.length.colonization + (1|Block)"
var <- "Root.diameter"

btaridity <- bootari(
  DF = dat,
  model = model_1,
  var = var,
  window = 115,
  NsimB = 500 )
colnames(btaridity) <- c("Intcp","TP","var")

model_11 <-"Total.biomass ~ Root.length.colonization+(1|Block)"
var1<-"Root.diameter"
btaridity1 <- bootari(
  DF = dat,
  model = model_11,
  var = var1,
  window = 115,
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

#write.csv(df1, "1_219_MW_115_TB_RLC_RAD.csv", quote=F)

df1=read.csv("1_219_MW_115_TB_RLC_RAD.csv",head=T,row.name=1)

df_mean <- df1 %>%
  dplyr::group_by(var) %>%
  dplyr::summarise(mean_TP = mean(TP, na.rm = TRUE))

fit_model <- lm(mean_TP ~ poly(var, 2, raw = TRUE), data = df_mean)
data.pred <- data.frame( mean_TP = predict(fit_model, newdata = data.frame(var = df1$var)), var = df1$var)

df1$group <- ifelse(
  df1$var == 0.260245689655172 ,
 "no", "sig" )

figs7e <- ggplot() + 
  geom_point(data = df1, aes(x = var, y = TP, color = group ), size = 2.5, alpha = 0.08, show.legend = F ) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "solid") +  
  labs(x = "Root diameter (mm)", y = "Standardized correlation effect") +
  theme(plot.background = element_blank(), panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() + 
  theme(axis.text = element_text(size = 11), axis.text.y = element_text(angle = 90, hjust = 0.5), axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12)) +
  theme(legend.position = c(0.55,0.95), legend.title = element_blank(), legend.box.background = element_blank(), legend.background = element_blank(),
        legend.text = element_text(size = 12), legend.key = element_blank()) +
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  geom_line(data = data.pred, aes(x = var, y = mean_TP, colour = "TB ~ Root length colonization"), size = 0.8, linetype = "solid",  show.legend = T) +  
  scale_color_manual(name = "", values = c("TB ~ Root length colonization" = "black","no" = "#386CB0", "sig" = "#ff8a00"), breaks = c("TB ~ Root length colonization"), labels = c("TB ~ Root length colonization" = "TB ~ Root length colonization")  )

figs7e


figs7g <- ggplot(df1, aes( x=as.factor(round(var, 3)),y= TP )) + 
  geom_boxplot(width=0.5,outlier.shape = 1,outlier.size = .3,)+ geom_hline(yintercept=0)+
  scale_color_manual(name = "", values = c("no" = "#386CB0", "sig" = "#ff8a00"))+
  geom_jitter(aes(color = group ), shape=1, size=.3,position=position_jitter(0.05),alpha=0.2, show.legend = F)+ 
  xlab("Root diameter (mm)")+ylab("")+ ggtitle("Standardized bootstrapped coefficients of correlation between TB and root length colonization") + theme_bw()+theme(panel.grid=element_blank())+theme(plot.background = element_blank())+
  theme(axis.text.x = element_text(size = 9,color = "black",angle = 60, hjust = 1),axis.ticks.x = element_blank(),title = element_text(size = 11)) +
  theme(axis.text.y = element_text(size = 9,color = "black"))+theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 11)) +
  scale_x_discrete(labels = function(x) format(as.numeric(x), nsmall = 3))+
  annotate("text", 
           x = as.factor(round(unique(df1$var[!(
             df1$var %in% c( 0.260245689655172
             )
           )]), 3)), 
           y = -0.32,  
           label = "*", 
           size = 5, 
           color = "black")

figs7g


pp <- cowplot::plot_grid( figs7a,  figs7b, figs7c,  figs7d, figs7e, 
                 ncol = 3, nrow = 2, labels = c('(a)', "(b)", "(c)",'(d)', "(e)"), 
                 label_fontface = "bold", label_x = -0.03, label_y = 1.02 , align = c("hv"))
pp

ppp1 <- cowplot::plot_grid( pp, legend, align = "", axis = "tpl", rel_widths = c( 1.6, 0.28 )) 
ppp1


ppp2 <- cowplot::plot_grid( figs7f,  figs7g,
                   ncol = 1, nrow = 2, labels = c( "(f)", "(g)"), 
                   label_fontface = "bold", label_x = -0.005, label_y = 1.02 , align = c("hv"))
ppp2

# 2_figs7_13*13
figs7 <- cowplot::plot_grid(
  ppp1, 
  ppp2, 
  ncol = 1,         
  align = "hv",      
  axis = "tb",      
  rel_heights = c(1, 0.7) 
)

figs7



######## Plant’s responses to mycorrhizal inoculation ----------------
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
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/greenhouse experiments")  
data = read.csv("greenhouse_data_2.csv",header=TRUE, row.names=1)

##
tree_raw <- read.tree("2022myco_45sp.nwk") 
plot(tree_raw)
rooted_tree1 <- root(tree_raw, outgroup = "Juncus", node = NULL, resolve.root = T, interactive = FALSE, edgelabel = FALSE)
plot(rooted_tree1)
rooted_tree  <- drop.tip(rooted_tree1, "Juncus")
plot(rooted_tree)

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

cols <- c("Amaryllidaceae"= "#6C8FA6", "Cyperaceae" = "#5F97C6", "Poaceae" = "#8FB9D1", 
          "Ranunculaceae" = "#BFA87F","Onagraceae" = "#B86A2E", "Brassicaceae"= "#B09AA8",  
          "Rosaceae" =  "#E0A16F", "Fabaceae" = "#FFA64D" , "Caryophyllaceae" = "#7FC1B5" , 
          "Polygonaceae"= "#6BA89F",  "Gentianaceae"= "#9e8c93","Lamiaceae"  = "#D94F30",
          "Plantaginaceae" = "#A83248", "Apiaceae" = "#B85A78", "Asteraceae" = "#E9746A") 


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

#1_fig3_10*9
fig3 <- p1 %>%
  insert_left(p00, width = 0.8) %>%          
  insert_right(blank, width = 0.03) %>%      
  insert_right(p2,   width = 1) %>%         
  insert_right(blank, width = 0.03) %>%     
  insert_right(p3,   width = 1)


fig3



######## Regression relationships of MPR and MGR with root length colonization, root diameter, and MRR ----------------
library(car)
library(ggplot2)
library(ggpmisc)
library(scales)
library(ggpubr)
library(cowplot)
library(ppcor)
library(doBy)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/greenhouse experiments")  
data = read.csv("greenhouse_data_2.csv",header=TRUE, row.names=1) 

cols <- c("Amaryllidaceae"= "#6C8FA6", "Cyperaceae" = "#5F97C6", "Poaceae" = "#8FB9D1", 
          "Ranunculaceae" = "#BFA87F","Onagraceae" = "#B86A2E", "Brassicaceae"= "#B09AA8",  
          "Rosaceae" =  "#E0A16F", "Fabaceae" = "#FFA64D" , "Caryophyllaceae" = "#7FC1B5" , 
          "Polygonaceae"= "#6BA89F",  "Gentianaceae"= "#9e8c93","Lamiaceae"  = "#D94F30",
          "Plantaginaceae" = "#A83248", "Apiaceae" = "#B85A78", "Asteraceae" = "#E9746A") 

se <- function(x) sd(x)/sqrt(length(x))
conf_95 <- function(x) t.test(x, conf.level = 0.95)$conf.int[2]

outlier_test = outlierTest(lm( MPR ~ Root.length.colonization , data = data)) 
outlier_test 
sumdata <- summaryBy( Root.length.colonization + MPR ~ Plant.species, data, FUN = c(mean, sd, se, conf_95))
sumdata
sumdata$Plant.family <- data$Plant.family[match(sumdata$Plant.species, data$Plant.species)]
sumdata$Root.length.colonization <- sumdata$Root.length.colonization.mean
sumdata$MPR <- sumdata$MPR.mean
colnames(sumdata)
fig4a <- ggplot(sumdata, aes(y = MPR , x = Root.length.colonization, color = Plant.family  ))+      
  scale_color_manual(values = cols)+  labs(color="Plant family")+
  scale_x_continuous(name ="Root length colonization (%)") + scale_y_continuous(name = "MPR" )+
  geom_errorbarh(aes(xmin = Root.length.colonization + Root.length.colonization.se, xmax = Root.length.colonization - Root.length.colonization.se), height=0, show.legend = F) +
  geom_errorbar(aes(ymin = MPR + MPR.se, ymax = MPR - MPR.se),  width=0, show.legend = F) +
  geom_point(size =  2.5, alpha = 1, shape = 16, aes(color = Plant.family), show.legend = F ) + 
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+  
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")), data = data, formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.98,label.y = 0.02, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  theme_bw() + theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), legend.text=element_text(size=11), legend.title = element_text(size=12, hjust=0.8),  axis.title=element_text(size=12)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))
fig4a


outlier_test = outlierTest(lm( MPR ~ Root.diameter , data = data)) 
outlier_test 
sumdata <- summaryBy( Root.diameter + MPR ~ Plant.species, data, FUN = c(mean, sd, se, conf_95))
sumdata
sumdata$Plant.family <- data$Plant.family[match(sumdata$Plant.species, data$Plant.species)]
sumdata$Root.diameter <- sumdata$Root.diameter.mean
sumdata$MPR <- sumdata$MPR.mean
colnames(sumdata)
fig4c <- ggplot( sumdata, aes(y = MPR , x = Root.diameter, color = Plant.family  ))+      
  scale_color_manual(values = cols)+  labs(color="Plant family")+
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "MPR" )+
  geom_errorbarh(aes(xmin = Root.diameter + Root.diameter.se, xmax = Root.diameter - Root.diameter.se), height=0, show.legend = F) +
  geom_errorbar(aes(ymin = MPR + MPR.se, ymax = MPR - MPR.se),  width=0, show.legend = F) +
  geom_point(size =  2.5, alpha = 1, shape = 16, aes(color = Plant.family), show.legend = F ) + 
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")), data = data, formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.02,label.y = 0.98, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  theme_bw() + theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), legend.text=element_text(size=11), legend.title = element_text(size=12, hjust=0.8),  axis.title=element_text(size=12)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))
fig4c


outlier_test = outlierTest(lm( MPR ~ MRR , data = data)) 
outlier_test 
sumdata <- summaryBy( MRR + MPR ~ Plant.species, data, FUN = c(mean, sd, se, conf_95))
sumdata
sumdata$Plant.family <- data$Plant.family[match(sumdata$Plant.species, data$Plant.species)]
sumdata$MRR <- sumdata$MRR.mean
sumdata$MPR <- sumdata$MPR.mean
colnames(sumdata)
fig4e <- ggplot(sumdata, aes(y = MPR , x = MRR, color = Plant.family  ))+      
  scale_color_manual(values = cols)+  labs(color="Plant family")+
  scale_x_continuous(name ="MRR") + scale_y_continuous(name = "MPR" )+
  geom_errorbarh(aes(xmin = MRR + MRR.se, xmax = MRR - MRR.se),  height=0, show.legend = F) +
  geom_errorbar(aes(ymin = MPR + MPR.se, ymax = MPR - MPR.se),  width=0, show.legend = F) +
  geom_point(size =  2.5, alpha = 1, shape = 16, aes(color = Plant.family), show.legend = F ) + 
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")), data = data, formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.98,label.y = 0.02, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  theme_bw() + theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), legend.text=element_text(size=11), legend.title = element_text(size=12, hjust=0.8),  axis.title=element_text(size=12)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))
fig4e


outlier_test = outlierTest(lm( MGR ~ Root.length.colonization , data = data)) 
outlier_test 
sumdata <- summaryBy( Root.length.colonization + MGR ~ Plant.species, data, FUN = c(mean, sd, se, conf_95))
sumdata
sumdata$Plant.family <- data$Plant.family[match(sumdata$Plant.species, data$Plant.species)]
sumdata$Root.length.colonization <- sumdata$Root.length.colonization.mean
sumdata$MGR <- sumdata$MGR.mean
colnames(sumdata)
fig4b <- ggplot(sumdata, aes(y = MGR , x = Root.length.colonization, color = Plant.family  ))+      
  scale_color_manual(values = cols)+  labs(color="Plant family")+
  scale_x_continuous(name ="Root length colonization (%)") + scale_y_continuous(name = "MGR" )+
  geom_errorbarh(aes(xmin = Root.length.colonization + Root.length.colonization.se, xmax = Root.length.colonization - Root.length.colonization.se), height=0, show.legend = F) +
  geom_errorbar(aes(ymin = MGR + MGR.se, ymax = MGR - MGR.se), width=0, show.legend = F ) +
  geom_point(size =  2.5, alpha = 1, shape = 16, aes(color = Plant.family), show.legend = F ) + 
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")), data = data, formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.02,label.y = 0.98, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  theme_bw() + theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), legend.text=element_text(size=11), legend.title = element_text(size=12, hjust=0.15),  
                     legend.position = "right", axis.title=element_text(size=12),legend.direction = "horizontal", legend.box = "horizontal") + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", nrow= 12, label.hjust = 0))
fig4b

outlier_test = outlierTest(lm( MGR ~ Root.diameter , data = data)) 
outlier_test 
sumdata <- summaryBy( Root.diameter + MGR ~ Plant.species, data, FUN = c(mean, sd, se, conf_95))
sumdata
sumdata$Plant.family <- data$Plant.family[match(sumdata$Plant.species, data$Plant.species)]
sumdata$Root.diameter <- sumdata$Root.diameter.mean
sumdata$MGR <- sumdata$MGR.mean
colnames(sumdata)

fig4d <- ggplot(sumdata, aes(y = MGR , x = Root.diameter, color = Plant.family  ))+      
  scale_color_manual(values = cols)+  labs(color="Plant family")+
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "MGR" )+
  geom_errorbarh(aes(xmin = Root.diameter + Root.diameter.se, xmax = Root.diameter - Root.diameter.se), height=0, show.legend = F) +
  geom_errorbar(aes(ymin = MGR + MGR.se, ymax = MGR - MGR.se), width=0, show.legend = F) +
  geom_point(size =  2.5, alpha = 1, shape = 16, aes(color = Plant.family), show.legend = F ) + 
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")), data = data, formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.02,label.y = 0.98, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  theme_bw() + theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), legend.text=element_text(size=11), legend.title = element_text(size=12, hjust=0.8),  axis.title=element_text(size=12)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", nrow= 12, label.hjust = 0))
fig4d


outlier_test = outlierTest(lm( MGR ~ MRR , data = data)) 
outlier_test 
sumdata <- summaryBy( MRR + MGR ~ Plant.species, data, FUN = c(mean, sd, se, conf_95))
sumdata
sumdata$Plant.family <- data$Plant.family[match(sumdata$Plant.species, data$Plant.species)]
sumdata$MRR <- sumdata$MRR.mean
sumdata$MGR <- sumdata$MGR.mean
colnames(sumdata)
fig4f <- ggplot( sumdata, aes( y = MGR , x = MRR, color = Plant.family  ))+      
  scale_color_manual(values = cols, breaks = c("Asteraceae",  "Apiaceae", "Plantaginaceae", 
                                               "Caryophyllaceae", "Polygonaceae", "Fabaceae", "Onagraceae", "Rosaceae",
                                               "Brassicaceae", "Ranunculaceae", "Poaceae", "Cyperaceae" ) )+  labs(color="Plant family")+
  scale_x_continuous(name ="MRR") + scale_y_continuous(name = "MGR" )+
  geom_errorbarh(aes(xmin = MRR + MRR.se, xmax = MRR - MRR.se), height =0, show.legend = F) +
  geom_errorbar(aes(ymin = MGR + MGR.se, ymax = MGR - MGR.se),  width=0, show.legend = F ) +
  geom_point(size =  2.5, alpha = 1, shape = 16, aes(color = Plant.family), show.legend = F ) + 
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")), data = data, formula = y ~ poly(x, 1, raw = TRUE), parse = TRUE,label.x = 0.98,label.y = 0.02, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  theme_bw() + theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), legend.text=element_text(size=11), legend.title = element_text(size=12, hjust=0.15),  axis.title=element_text(size=12)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", nrow= 12, label.hjust = 0))
fig4f


legend <- get_legend(fig4f) 


pp1 <- plot_grid( fig4a, fig4c, fig4e,  fig4b,  fig4d,  fig4f, 
                  ncol = 3, nrow = 2, labels = c('(a)', "(c)", "(e)",'(b)', "(d)", "(f)"), 
                  label_fontface = "bold", label_x = -0.03, label_y = 1.02 , align = c("hv"))
pp1


#1_fig4_12*8
fig4 <- cowplot::plot_grid( pp1, legend, align = "", axis = "tpl", rel_widths = c( 1.7, 0.28 )) 
fig4


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


######## Phylogenetic generalized least squares regressions ----------------
library(car)
library(caper) 
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/greenhouse experiments")  
data = read.csv("greenhouse_data_mean.csv",header=TRUE, row.names=1) 

tree_raw <- read.tree("2022myco_45sp.nwk") 
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

data_s = data[-c(11,12),]
outlier_test = outlierTest(lm(P.concentrations ~ Root.length.colonization, data =  data_s)) 
outlier_test
comp_data_s <- comparative.data( phy = rooted_tree, data = data_s , names.col = "Plant.species", vcv = TRUE, vcv.dim = 2)
pgls_model <- pgls(P.concentrations ~ Root.length.colonization, data = comp_data_s, lambda = "ML") 
summary(pgls_model)

data_s = data[-c(30,32),]
outlier_test = outlierTest(lm(Total.biomass ~ Root.length.colonization, data =  data_s)) 
outlier_test
comp_data_s <- comparative.data( phy = rooted_tree, data = data_s , names.col = "Plant.species", vcv = TRUE, vcv.dim = 2)
pgls_model <- pgls(Total.biomass ~ Root.length.colonization, data = comp_data_s, lambda = "ML") 
summary(pgls_model)

###
outlier_test = outlierTest(lm( MPR ~ Root.length.colonization, data =  data)) 
outlier_test
pgls_model <- pgls( MPR ~ Root.length.colonization , data = comp_data, lambda = "ML") 
summary(pgls_model)

outlier_test = outlierTest(lm( MPR ~ Root.diameter, data =  data)) 
outlier_test
pgls_model <- pgls( MPR ~ Root.diameter , data = comp_data, lambda = "ML")
summary(pgls_model)

outlier_test = outlierTest(lm( MPR ~ MRR, data =  data)) 
outlier_test
pgls_model <- pgls( MPR ~ MRR , data = comp_data, lambda = "ML") 
summary(pgls_model)


###
outlier_test = outlierTest(lm( MGR ~ Root.length.colonization, data =  data)) 
outlier_test
pgls_model <- pgls( MGR ~ Root.length.colonization , data = comp_data, lambda = "ML") 
summary(pgls_model)

outlier_test = outlierTest(lm( MGR ~ Root.diameter, data =  data)) 
outlier_test
pgls_model <- pgls( MGR ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)

outlier_test = outlierTest(lm( MGR ~ MRR, data =  data)) 
outlier_test
pgls_model <- pgls( MGR ~ MRR , data = comp_data, lambda = "ML") 
summary(pgls_model)

