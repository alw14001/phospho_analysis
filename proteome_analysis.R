---
  title: "Proteome Analysis"
author: "Timothy E. Moore"
date: "`r Sys.Date()`"
output:
  rmdformats::readthedown:
  highlight: kate
code_folding: hide
editor_options: 
  chunk_output_type: console
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

##msDiaLogue- PROCESSING STEPS###################################

```{r}

#load packages
library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)
library(msDiaLogue)
library(PCAtest)
library(FactoMineR)


#load proteome dataset and run preprocessing step to prepare data
fileName <- "audrey_rawdata.csv"
dataSet <- preprocessing(fileName)

#perform log based transformation to stabilize variance
dataTrans <- transform(dataSet, logFold = 2)

#prints name column from dataTrans
#names(dataTrans)

#plot PCA
res.pca <- PCA(dataTrans[, c(1,4:9410)], 
               quali.sup = 1, graph = TRUE)
dataTrans$pc1 <-res.pca$ind$coord[, 1]
dataTrans$pc2 <-res.pca$ind$coord[, 2]

#plot a nicer version of PCA
ggplot(data = dataTrans, aes(x = pc1, y = pc2, color = R.Condition, 
                             label = R.Replicate))+
  #  geom_point()+
  geom_text()+
  geom_vline(xintercept = 0, lty =2)+
  geom_hline(yintercept = 0, lty = 2)


#normalize data with boxplots showing after normalization
dataNorm <- normalize(dataTrans, normalizeType = "mean")

#imputation using the sequential knn algorithm (kim et al 2024)
dataImput <- impute.knn_seq(dataNorm)

#remove any missing values that remain after imputation
dataImput <- filterNA(dataImput, saveRm = TRUE)

#computes summary statistics for each condition for every protein
dataSumm <- summarize(dataImput, saveSumm = TRUE)

```
##ANALYSIS#######################################################

### Do all pairwise comparisons

```{r}

pairwise_comparisons <- combn(unique(dataImput$R.Condition), 
                              2, simplify = TRUE) %>%
  t() %>% as.data.frame() %>% arrange(desc(V1))

pairwise_analysis <- list()
pairwise_visualize <-list()

for (i in 1: nrow(pairwise_comparisons)){
  
  pairwise_analysis[[i]] <- analyze(dataImput, testType = "mod.t-test", 
                                    cond = pairwise_comparisons[i, ])
  
  pairwise_visualize[[i]] <- visualize(pairwise_analysis[[i]], 
                                       graphType = "volcano", 
                                       P.thres = 0.05, logF.thres = 0.5)+
    ggtitle(paste(pairwise_comparisons[i, 1], pairwise_comparisons[i, 2],
                  sep = " & "))
  print(i)  
}

pairwise_visualize[1]

P.thres = 0.05
logF.thres = 0.5


p7.a <- analyze(dataImput, testType = "mod.t-test", 
                cond = c("p7het", "p7wt"))
p14.a <- analyze(dataImput, testType = "mod.t-test", 
                 cond = c("p14het", "p14wt"))
p28.a <- analyze(dataImput, testType = "mod.t-test", 
                 cond = c("p28het", "p28wt"))
p7p14.a <- analyze(dataImput, testType = "mod.t-test", 
                   cond = c("p7het", "p14wt"))


plotData <- data.frame(t(p28.a)) %>%
  mutate(Significant = case_when(
    P.value < P.thres & Difference > logF.thres ~ "Up",
    P.value < P.thres & Difference < -logF.thres ~ "Down",
    P.value < P.thres & Difference >= -logF.thres & Difference <= logF.thres ~ "Inconclusive",
    P.value >= P.thres ~ "No",
    TRUE ~ "Unknown")) %>% ## optional catch-all for other cases
  mutate(delabel = ifelse(! Significant %in% c("No", "Inconclusive"),
                          gsub("_.*", "", colnames(dataSet)), NA))

ggplot(plotData, aes(x = Difference, y = -log10(P.value),
                     col = Significant)) +
  geom_vline(xintercept = c(-logF.thres, logF.thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(P.thres), linetype = "dashed") +
  geom_point() +
  xlim(-3.5, 3.5) +
  scale_color_manual(values = c("Down" = "blue", "Up" = "red",
                                "Inconclusive" = "gray", "No" = "gray20")) +
  labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"P-value")) +
  theme_classic() +
  theme(legend.position = c(.9,.8), plot.title = element_text(hjust = 0.5))




p7.v <- visualize(p7.a,
                  graphType = "volcano", 
                  P.thres = 0.05, logF.thres = 0.5)
p14.v <- visualize(p14.a,
                   graphType = "volcano", 
                   P.thres = 0.05, logF.thres = 0.5)
p28.v <- visualize(p28.a,
                   graphType = "volcano", 
                   P.thres = 0.05, logF.thres = 0.5)
p7p14.v <- visualize(p7p14.a,
                     graphType = "volcano", 
                     P.thres = 0.05, logF.thres = 0.5)


p7.v+ggtitle("p7")+geom_text_repel(size=1)
p14.v+ggtitle("p14")  
p28.v+ggtitle("p28")  
p7p14.v+ggtitle()



```

## volcano plots

```{r}
library(patchwork)
# pairwise_comparisons
patchwork::wrap_plots(p7.v, p14.v, p28.v)
ggsave("volcano_plot.png", height = 8, width = 12, units = "in", 
       dpi = 400)

p7p14.v
ggsave("volcano_plot_p7p14.png", height = 6, width = 7, units = "in", 
       dpi = 400)


dev.new()
patchwork::wrap_plots(pairwise_visualize[c(1, 4, 10, 
                                           13,5, 14)], 
                      nrow = 2, ncol = 3)



dataImput$R.Condition
dataNorm$A12345_MOUSE.B7ZBV9_MOUSE

dataNorm2 <- dataNorm %>% arrange(R.Condition)
dataNorm2$R.Condition
dataNorm2$time <- c(rep(14, 6), rep(28, 6), rep(7, 6))
dataNorm2$cond <- c(rep(c("het", "het", "het", "wt", "wt", "wt"), 3))

dataNorm2$R.Condition, A12345_MOUSE.B7ZBV9_MOUSE

new_df <- dataNorm[,c("R.Condition","R.Replicate", "A12345_MOUSE.B7ZBV9_MOUSE", "KCNQ3_MOUSE", "E9Q9F1_MOUSE.G3UYG5_MOUSE.KCNQ5_MOUSE")]
new_df

colnames(new_df) <- c('Condition', 'Replicate', 'Kcnq2', 'Kcnq3', 'Kcnq5')
new_df
df2 %>% dplyr::select(dataNorm$R.Condition, dataNorm$R.Replicate, dataNorm$A12345_MOUSE.B7ZBV9_MOUSE)

write.csv(new_df, '~/Documents/GitHub/MsDialogue/kcnq_values.csv')

write.csv(dataImput,  '~/Documents/GitHub/MsDialogue/new_normalized_data.csv')
dataNorm2[, c("R.Condition", "time", "cond")]

dat.sum <- dataNorm2 %>% group_by(cond, time) %>%
  summarise(mean.p =  mean(KCNQ3_MOUSE, na.rm = TRUE), 
            sd.p = sd(KCNQ3_MOUSE, na.rm = TRUE))

ggplot(data = dat.sum, aes(x = time, y = mean.p, color = cond,group=cond))+
  geom_point(alpha = 0.2)+
  geom_point(data = dat.sum, aes(x = time, y = mean.p, group = cond, color = cond), size = 3)+
  geom_errorbar(data = dat.sum, aes(ymin = mean.p - sd.p, ymax = mean.p + sd.p, group = cond), width = 0.3)+
  theme_classic()+
  scale_color_manual(values = c("darkblue", "darkred"))+
  xlab("Age (Days)")+
  ylab("KCNQ3 Expression")
ggsave("lineplot.png", height = 6, width = 6.5, units = "in", 
       dpi= 400)




```

#---------------------------------------------#
# PCA Analysis

```{r}

library(FactoMineR)
library(factoextra)
library(dplyr)
library(ggplot2)

dat.o <- read.csv("audrey_rawdata.csv", header = TRUE)
clust <- read.csv("Protein_clusters.csv", header = TRUE)
head(clust)

clust %>% dplyr::select(PG.ProteinAccessions, contains("KCNQ2"))

dat.clust <- dat.o[dat.o$PG.Genes %in% clust$PG.Genes,]

colnames(dataImput)[colnames(dataImput) %in% dat.clust$PG.ProteinNames]

```

PG.Protein Names NOT present in the imputation data

```{r}
# PG.Protein Names NOT present in the imputation data
unique(dat.clust$PG.ProteinNames)[!unique(dat.clust$PG.ProteinNames) %in% colnames(dataImput)]
# [1] "A0A067XG51_MOUSE"            # single hit                            
# [2] "F7BA29_MOUSE.NMD3B_MOUSE"    # single hit                     
# [3] "S4R1M4_MOUSE"                # single hit                          
# [4] "A0A087WRP9_MOUSE" also Ank3? # removed at imputation               
# [5] "S4R1S2_MOUSE" also Ank3?     # removed at imputation                 
# [6] "A0A0J9YVF2_MOUSE"            # single hit                      
# [7] "A2A542_MOUSE"                # some single hits                      
# [8] "E0CX94_MOUSE"  also Kcnip2?  # some single hits                      
# [9] "E0CXC0_MOUSE"  also Kcnip2?  # some single hits                      
# [10] "A0A1D5RLU3_MOUSE"           # some single hits                       
# [11] "A0A286YCQ2_MOUSE.A0A286YDR7_MOUSE.A0A286YDU9_MOUSE" # some single hits
# [12] "A0A571BFX6_MOUSE"           # single hits

```


```{r}
# get protein names of proteins associated with K genes
prot.k <- unique(dat.clust$PG.ProteinNames[dat.clust$PG.ProteinAccessions %in% 
                                             clust[clust$Channel.group=="K+ Channels", ]$PG.ProteinAccessions])

prot.ca <- unique(dat.clust$PG.ProteinNames[dat.clust$PG.ProteinAccessions %in% 
                                              clust[clust$Channel.group=="Ca2+ Channels", ]$PG.ProteinAccessions])
prot.nmda <- unique(dat.clust$PG.ProteinNames[dat.clust$PG.ProteinAccessions %in%                 clust[clust$Channel.group=="NMDA/AMPA",]$PG.ProteinAccessions])
prot.hcn <- unique(dat.clust$PG.ProteinNames[dat.clust$PG.ProteinAccessions %in% 
                                               clust[clust$Channel.group=="Ankyrin/HCN", ]$PG.ProteinAccessions])


# filter imputed data to those specific proteins
focal.K <- dataImput %>% dplyr::select(R.Condition, 
                                       all_of(prot.k[prot.k%in%colnames(dataImput)]))
focal.CA <- dataImput %>% dplyr::select(R.Condition, 
                                        all_of(prot.ca[prot.ca%in%colnames(dataImput)]))
focal.NMDA <- dataImput %>% dplyr::select(R.Condition,                                                all_of(prot.nmda[prot.nmda%in%colnames(dataImput)]))
focal.HCN <- dataImput %>% dplyr::select(R.Condition, 
                                         all_of(prot.hcn[prot.hcn%in%colnames(dataImput)]))
```


```{r}
# run obs data PCA

names(dataImput)


res.pca.K <- PCA(focal.K, 
                 quali.sup = 1, graph = TRUE)
res.pca.K <- PCA(focal.K, 
                 quali.sup = 1, graph = TRUE)

res.pca.K$var$coord
res.pca.K$eig[1, 1]/sum(res.pca.K$eig[, 1])

focal.K$pc1 <- res.pca.K$ind$coord[, 1]
focal.K$pc2 <- res.pca.K$ind$coord[, 2]

ggplot(data = focal.K, aes(x = pc1, y = pc2, color =R.Condition))+
  geom_point()+
  geom_vline(xintercept = 0, lty=2)+
  geom_hline(yintercept = 0, lty =2)


res.pca.CA <- PCA(focal.CA, 
                  quali.sup = 1, graph = TRUE)

res.pca.NMDA <- PCA(focal.NMDA, 
                    quali.sup = 1, graph = TRUE)

res.pca.HCN <- PCA(focal.HCN, 
                   quali.sup = 1, graph = FALSE)



# calculate observed distances

obs.co.K <-res.pca.K$quali.sup$coord[, 1:2] %>% as.data.frame() %>% 
  tibble::rownames_to_column("cond") %>% mutate(group="K")
obs.co.CA <-res.pca.CA$quali.sup$coord[, 1:2] %>% as.data.frame() %>% 
  tibble::rownames_to_column("cond") %>% mutate(group="CA")
obs.co.NMDA <-res.pca.NMDA$quali.sup$coord[, 1:2] %>% as.data.frame() %>% 
  tibble::rownames_to_column("cond") %>% mutate(group="NMDA")
obs.co.HCN <-res.pca.HCN$quali.sup$coord[, 1:2] %>% as.data.frame() %>% 
  tibble::rownames_to_column("cond") %>% mutate(group="HCN")

obs.co.comb <- rbind(obs.co.K, obs.co.CA, obs.co.NMDA, obs.co.HCN)

combos <- combn(sort(obs.co.K$cond, decreasing=TRUE), 2, simplify = TRUE) %>% t() %>%   
  data.frame() 

combos[6, ]<-combos[6, c(2, 1)]
combos[8, ]<-combos[8, c(2, 1)]
combos[13, ]<-combos[13, c(2, 1)]

combos$distance <- NA

for (i in 1:15){
  i=1
  x1 <- obs.co.K[obs.co.K$cond==combos[i, "X1"], "Dim.1"]
  x2 <- obs.co.K[obs.co.K$cond==combos[i, "X1"], "Dim.1"]
  y1 <- obs.co.K[obs.co.K$cond==combos[i, "X1"], "Dim.2"]
  y2 <- obs.co.K[obs.co.K$cond==combos[i, "X1"], "Dim.2"]
  dist.obs <- sqrt((x1 - x2)^2 + (y1 - y2)^2)
  
}


```


```{r}
# generate null PCAs and distances
null.dist.mats.K <- list()
null.dist.mats.CA <- list()
null.dist.mats.NMDA <- list()
null.dist.mats.HCN <- list()

null.pca.K <- list()
null.pca.CA <- list()
null.pca.NMDA <- list()
null.pca.HCN <- list()

for(i in 1:1000){
  
  nulls.K <- sample(names(dataImput)[4:length(names(dataImput))], 
                    ncol(focal.K)-1)
  null.pca.K <- PCA(dataImput %>% dplyr::select(R.Condition, all_of(nulls.K)), 
                    quali.sup = 1, graph = FALSE)
  null.dist.mats.K[[i]] <- dist(null.pca.K$quali.sup$coord[, 1:2])
  
  nulls.CA <- sample(names(dataImput)[4:length(names(dataImput))], 
                     ncol(focal.CA)-1)
  null.pca.CA <- PCA(dataImput %>% dplyr::select(R.Condition, all_of(nulls.CA)), 
                     quali.sup = 1, graph = FALSE)
  null.dist.mats.CA[[i]] <- dist(null.pca.CA$quali.sup$coord[, 1:2])
  
  nulls.NMDA <- sample(names(dataImput)[4:length(names(dataImput))], 
                       ncol(focal.NMDA)-1)
  null.pca.NMDA <- PCA(dataImput %>% dplyr::select(R.Condition, 
                                                   all_of(nulls.NMDA)), 
                       quali.sup = 1, graph = FALSE)
  null.dist.mats.NMDA[[i]] <- dist(null.pca.NMDA$quali.sup$coord[, 1:2])
  
  nulls.HCN <- sample(names(dataImput)[4:length(names(dataImput))], 
                      ncol(focal.HCN)-1)
  null.pca.HCN <- PCA(dataImput %>% dplyr::select(R.Condition, all_of(nulls.HCN)), 
                      quali.sup = 1, graph = FALSE)
  null.dist.mats.HCN[[i]] <- dist(null.pca.HCN$quali.sup$coord[, 1:2])
  
}

null.dist.m.K <- list()
null.dist.m.CA <- list()
null.dist.m.NMDA <- list()
null.dist.m.HCN <- list()

for(i in 1:1000){
  m.K <- null.dist.mats.K[[i]] %>% as.matrix()
  m.CA <- null.dist.mats.CA[[i]] %>% as.matrix()
  m.NMDA <- null.dist.mats.NMDA[[i]] %>% as.matrix()
  m.HCN <- null.dist.mats.HCN[[i]] %>% as.matrix()
  
  m2.K <- reshape::melt(m.K)[reshape::melt(upper.tri(m.K))$value,]
  m2.CA <- reshape::melt(m.CA)[reshape::melt(upper.tri(m.CA))$value,]
  m2.NMDA <- reshape::melt(m.NMDA)[reshape::melt(upper.tri(m.NMDA))$value,]
  m2.HCN <- reshape::melt(m.HCN)[reshape::melt(upper.tri(m.HCN))$value,]
  
  names(m2.K) <- c("cond1", "cond2", "distance")
  names(m2.CA) <- c("cond1", "cond2", "distance")
  names(m2.NMDA) <- c("cond1", "cond2", "distance")
  names(m2.HCN) <- c("cond1", "cond2", "distance")
  
  m2.K$rep <- i
  m2.CA$rep <- i
  m2.NMDA$rep <- i
  m2.HCN$rep <- i
  
  null.dist.m.K[[i]] <- m2.K
  null.dist.m.CA[[i]] <- m2.CA
  null.dist.m.NMDA[[i]] <- m2.NMDA
  null.dist.m.HCN[[i]] <- m2.HCN
}

#merge all simulated null distances into data frame
null.dist.df.K <- do.call("rbind", null.dist.m.K)
null.dist.df.CA <- do.call("rbind", null.dist.m.CA)
null.dist.df.NMDA <- do.call("rbind", null.dist.m.NMDA)
null.dist.df.HCN <- do.call("rbind", null.dist.m.HCN)

```


## plot results
see:
  https://cran.r-project.org/web/packages/ggdist/vignettes/slabinterval.html

```{r}
library(ggdist)
library(colorspace)
null.dist.df.K$comb <- paste(null.dist.df.K$cond1, null.dist.df.K$cond2,
                             sep = " & ")
dist.obs.df.K$comb <- paste(dist.obs.df.K$cond1, dist.obs.df.K$cond2,
                            sep = " & ")
null.dist.df.CA$comb <- paste(null.dist.df.CA$cond1, null.dist.df.CA$cond2,
                              sep = " & ")
dist.obs.df.CA$comb <- paste(dist.obs.df.CA$cond1, dist.obs.df.CA$cond2,
                             sep = " & ")

null.dist.df.NMDA$comb <- paste(null.dist.df.NMDA$cond1, null.dist.df.NMDA$cond2,
                                sep = " & ")
dist.obs.df.NMDA$comb <- paste(dist.obs.df.NMDA$cond1, dist.obs.df.NMDA$cond2,
                               sep = " & ")

null.dist.df.HCN$comb <- paste(null.dist.df.HCN$cond1, null.dist.df.HCN$cond2,
                               sep = " & ")
dist.obs.df.HCN$comb <- paste(dist.obs.df.HCN$cond1, dist.obs.df.HCN$cond2,
                              sep = " & ")


dist.obs.df.K <- dist.obs.df.K %>% arrange(distance)
dist.obs.df.CA <- dist.obs.df.CA %>% arrange(distance)
dist.obs.df.NMDA <- dist.obs.df.NMDA %>% arrange(distance)
dist.obs.df.HCN <- dist.obs.df.HCN %>% arrange(distance)


dist.obs.df.K$comb <- factor(dist.obs.df.K$comb, 
                             levels = dist.obs.df.K$comb, 
                             ordered = TRUE)
dist.obs.df.CA$comb <- factor(dist.obs.df.CA$comb, 
                              levels = dist.obs.df.CA$comb, 
                              ordered = TRUE)
dist.obs.df.NMDA$comb <- factor(dist.obs.df.NMDA$comb, 
                                levels = dist.obs.df.NMDA$comb, 
                                ordered = TRUE)
dist.obs.df.HCN$comb <- factor(dist.obs.df.HCN$comb, 
                               levels = dist.obs.df.HCN$comb, 
                               ordered = TRUE)

null.dist.df.K$comb <- factor(null.dist.df.K$comb, 
                              levels = dist.obs.df.K$comb, 
                              ordered = TRUE)
null.dist.df.CA$comb <- factor(null.dist.df.CA$comb, 
                               levels = dist.obs.df.CA$comb, 
                               ordered = TRUE)
null.dist.df.NMDA$comb <- factor(null.dist.df.NMDA$comb, 
                                 levels = dist.obs.df.NMDA$comb, 
                                 ordered = TRUE)
null.dist.df.HCN$comb <- factor(null.dist.df.HCN$comb, 
                                levels = dist.obs.df.HCN$comb, 
                                ordered = TRUE)


unique(null.dist.df.K$comb)
within <- c("p7het & p7wt", "p14het & p14wt", "p28het & p28wt")

ggplot(data = null.dist.df.K , 
       aes(y = distance, 
           x = comb))+
  stat_halfeye()+
  geom_point(data = dist.obs.df.K , 
             aes(y = distance,
                 x = comb), color = "darkred", 
             size = 4.5)+
  coord_flip()+
  theme_minimal()

#d  = sqrt((x1-x2)^2*wx + (y1-y2)^2*wy)

ggplot(data = null.dist.df.CA , aes(y = distance, 
                                    x = comb))+
  stat_halfeye()+
  geom_point(data = dist.obs.df.CA  , 
             aes(y = distance,
                 x = comb), color = "darkred", 
             size = 4.5)+
  coord_flip()+
  theme_minimal()


ggplot(data = null.dist.df.NMDA , aes(y = distance, 
                                      x = comb))+
  stat_halfeye()+
  geom_point(data = dist.obs.df.NMDA  , 
             aes(y = distance,
                 x = comb), color = "darkred", 
             size = 4.5)+
  coord_flip()+
  theme_minimal()



```


##Compare differences (in distances) between conditions and time points for observed vs null samples

```{r}

dist.obs.df.K # data frame of observed distances
p14_dif_p7_dif <-dist.obs.df.K[1, 3] - dist.obs.df.K[2, 3]
dif_obs <-dist.obs.df.K[11, 3] - dist.obs.df.K[8, 3]


# data frame of null distances
test <- null.dist.df.K[(null.dist.df.K$cond1=="p14wt" & null.dist.df.K$cond2=="p7wt")| 
                         (null.dist.df.K$cond1=="p14het" & null.dist.df.K$cond2=="p7het"),]

test.l <- test %>% dplyr::select(distance, comb, rep) %>% pivot_wider(values_from = distance,
                                                                      names_from = comb, 
                                                                      id_cols = rep) %>% 
  mutate(difference = `p14wt & p7wt` - `p14het & p7het`)

head(test.l)


ggplot()+
  geom_histogram(data = test.l, aes(x = difference))+
  geom_vline(xintercept = dif_obs, lty = 3, col = "red")

```


g <- visualize(dat.sum, type = "heatmap")
dat.sum <- dataNorm2 %>% group_by(R.Condition, vars) %>%
  summarise_at(vars("KCNQ3_MOUSE", "A12345_MOUSE.B7ZBV9_MOUSE"), mean) 

dataNorm2 %>%
  group_by(R.Condition) %>%
  summarise_at(vars(-Month), funs(mean(., na.rm=TRUE)))

dat.sum

