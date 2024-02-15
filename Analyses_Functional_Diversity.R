#1. Functional space and hypervolume estimation by region
library(pacman)
pacman::p_load(dplyr, plyr, readr, tbible, FD, ade4, cowplot, mice, reshape2, tidyr, ks, hypervolume, alphallhu, purrr, TTR, plotrix, agricolae, psych)

library(factoextra)
library(ggrepel)
library(tibble)

#### Load data

# load: trait data
trait <- read.csv("C:\\Users\\camil\\OneDrive - SELVA Investigación para la conservación en el neotropico\\Documentos\\Camila\\Publicaciones\\Chapman_Historical Occupancy\\All_species_traits_CorrJuli_FINAL.csv", sep = ";")

head(trait)
# Exploring trait correlations (verify trait column numbers)

just_trait <- trait[, c(48,49,50,51,52,53,54,55,56,57,58)]
head(just_trait)
res <- cor(just_trait)
round(res, 2)

library(corrplot)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

library("PerformanceAnalytics")
chart.Correlation(just_trait, histogram=TRUE, pch=19)

#Select combination of traits with lower correlations


#Single PCA for EACH REGION

#z Transform traits
# function to z-transform data
scale_z <- function(x){
  (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
}

### add scaled trait measures to traits table
trait_z <- trait %>% 
  mutate(log_weight = log(Mass)) %>% 
  mutate(log_Beak.Length_Culmen = log(Beak.Length_Culmen)) %>% 
  mutate(log_Beak.Width = log(Beak.Width)) %>% 
  mutate(log_Wing.Length = log(Wing.Length)) %>% 
  mutate(log_Tail.Length = log(Tail.Length)) %>% 
  mutate(log_Tarsus.Length = log(Tarsus.Length)) %>% 
  mutate(weight_z = scale_z(log_weight)) %>% 
  mutate(T_cul_z = scale_z(log_Beak.Length_Culmen)) %>% 
  mutate(Bill_w_z = scale_z(log_Beak.Width)) %>% 
  mutate(Wing_l_z = scale_z(log_Wing.Length)) %>% 
  mutate(Tail_l_z = scale_z(log_Tail.Length)) %>% 
  mutate(Tars_z = scale_z(log_Tarsus.Length)) %>% 
  mutate(HWI_z = scale_z(Hand.Wing.Index)) %>% 
  as.data.frame

#subset scaled measures
trait_z <- trait_z %>% 
  select(Species,MODERN, HISTORICAL,MIR_H_corr,MIR_M,MIR_EXT,MIR_NEW, BARB_H_corr,BARB_M,BARB_EXT,BARB_NEW,SAGU_H_corr, SAGU_M,SAGU_EXT,SAGU_NEW, TOCHE_H_corr,TOCHE_M,TOCHE_EXT, TOCHE_NEW,HON_H_corr,HON_M,HON_EXT, HON_NEW,FLOR_H_corr, FLOR_M,FLOR_EXT, FLOR_NEW, FUSA_H_corr, FUSA_M,FUSA_EXT, FUSA_NEW,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) 

MIR <- trait_z %>% 
  filter(MIR_H_corr >= 0)

BARB <- trait_z %>% 
  filter(BARB_H_corr >= 0)

SAGU <- trait_z %>% 
  filter(SAGU_H_corr >= 0)

TOCHE <- trait_z %>% 
  filter(TOCHE_H_corr >= 0)

HON <- trait_z %>% 
  filter(HON_H_corr >= 0)

FLOR <- trait_z %>% 
  filter(FLOR_H_corr >= 0)

FUSA <- trait_z %>% 
  filter(FUSA_H_corr >= 0)
#### --------------------------------------------------------------
## PCA ## MIRAFLORES
#### --------------------------------------------------------------

#table just with columns of morphological traits and species names as rownames
tr_mi_z_MIR <- MIR %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)%>% 
  column_to_rownames(var = "Species")

#Run PCA
pcaTotal_MIR <- princomp(tr_mi_z_MIR, cor = TRUE, scores = TRUE)
summary(pcaTotal_MIR)

#Scores
scoresPCATotal_MIR <- as.data.frame(pcaTotal_MIR$scores) %>% 
  tibble::rownames_to_column("Species")

scoresPCATotal_MIR <- scoresPCATotal_MIR %>% 
  # convert long to wide
  tidyr::gather(key, value, -Species) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# loadings
loadingsPCATotal_MIR <- as.data.frame(unclass(pcaTotal_MIR$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()

loadingsPCATotal_MIR <- loadingsPCATotal_MIR %>% 
  # convert long to wide
  tidyr::gather(key, value, -trait) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# scalar to adjust arrow size
sc_mi <- 7

loadings_sc_MIR <- loadingsPCATotal_MIR %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")), funs(.*sc_mi)) %>% 
  # variable names
  mutate(trait = c("Bill width","Hand-Wing Index", "Total culmen", "Tail length", "Tarsus", "Body mass", "Wing length"))

#######################

#### --------------------------------------------------------------
## PCA ## BARBACOAS
#### --------------------------------------------------------------

#table just with columns of morphological traits and species names as rownames
tr_mi_z_BARB <- BARB %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)%>% 
  column_to_rownames(var = "Species")

#Run PCA
pcaTotal_BARB <- princomp(tr_mi_z_BARB, cor = TRUE, scores = TRUE)
summary(pcaTotal_BARB)

#Scores
scoresPCATotal_BARB <- as.data.frame(pcaTotal_BARB$scores) %>% 
  tibble::rownames_to_column("Species")

scoresPCATotal_BARB <- scoresPCATotal_BARB %>% 
  # convert long to wide
  tidyr::gather(key, value, -Species) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# loadings
loadingsPCATotal_BARB <- as.data.frame(unclass(pcaTotal_BARB$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()

loadingsPCATotal_BARB <- loadingsPCATotal_BARB %>% 
  # convert long to wide
  tidyr::gather(key, value, -trait) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# scalar to adjust arrow size
sc_mi <- 7

loadings_sc_BARB <- loadingsPCATotal_BARB %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")), funs(.*sc_mi)) %>% 
  # variable names
  mutate(trait = c("Bill width","Hand-Wing Index", "Total culmen", "Tail length", "Tarsus", "Body mass", "Wing length"))

#######################

#### --------------------------------------------------------------
## PCA ## SAN AGUSTIN
#### --------------------------------------------------------------

#table just with columns of morphological traits and species names as rownames
tr_mi_z_SAGU <- SAGU %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)%>% 
  column_to_rownames(var = "Species")

#Run PCA
pcaTotal_SAGU <- princomp(tr_mi_z_SAGU, cor = TRUE, scores = TRUE)
summary(pcaTotal_SAGU)

#Scores
scoresPCATotal_SAGU <- as.data.frame(pcaTotal_SAGU$scores) %>% 
  tibble::rownames_to_column("Species")

scoresPCATotal_SAGU <- scoresPCATotal_SAGU %>% 
  # convert long to wide
  tidyr::gather(key, value, -Species) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# loadings
loadingsPCATotal_SAGU <- as.data.frame(unclass(pcaTotal_SAGU$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()

loadingsPCATotal_SAGU <- loadingsPCATotal_SAGU %>% 
  # convert long to wide
  tidyr::gather(key, value, -trait) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# scalar to adjust arrow size
sc_mi <- 7

loadings_sc_SAGU <- loadingsPCATotal_SAGU %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")), funs(.*sc_mi)) %>% 
  # variable names
  mutate(trait = c("Bill width","Hand-Wing Index", "Total culmen", "Tail length", "Tarsus", "Body mass", "Wing length"))

#######################

#### --------------------------------------------------------------
## PCA ## TOCHE
#### --------------------------------------------------------------

#table just with columns of morphological traits and species names as rownames
tr_mi_z_TOCHE <- TOCHE %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)%>% 
  column_to_rownames(var = "Species")

#Run PCA
pcaTotal_TOCHE <- princomp(tr_mi_z_TOCHE, cor = TRUE, scores = TRUE)
summary(pcaTotal_TOCHE)

#Scores
scoresPCATotal_TOCHE <- as.data.frame(pcaTotal_TOCHE$scores) %>% 
  tibble::rownames_to_column("Species")

scoresPCATotal_TOCHE <- scoresPCATotal_TOCHE %>% 
  # convert long to wide
  tidyr::gather(key, value, -Species) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# loadings
loadingsPCATotal_TOCHE <- as.data.frame(unclass(pcaTotal_TOCHE$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()

loadingsPCATotal_TOCHE <- loadingsPCATotal_TOCHE %>% 
  # convert long to wide
  tidyr::gather(key, value, -trait) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# scalar to adjust arrow size
sc_mi <- 7

loadings_sc_TOCHE <- loadingsPCATotal_TOCHE %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")), funs(.*sc_mi)) %>% 
  # variable names
  mutate(trait = c("Bill width","Hand-Wing Index", "Total culmen", "Tail length", "Tarsus", "Body mass", "Wing length"))

#######################

#### --------------------------------------------------------------
## PCA ## HONDA
#### --------------------------------------------------------------

#table just with columns of morphological traits and species names as rownames
tr_mi_z_HON <- HON %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)%>% 
  column_to_rownames(var = "Species")

#Run PCA
pcaTotal_HON <- princomp(tr_mi_z_HON, cor = TRUE, scores = TRUE)
summary(pcaTotal_HON)

#Scores
scoresPCATotal_HON <- as.data.frame(pcaTotal_HON$scores) %>% 
  tibble::rownames_to_column("Species")

scoresPCATotal_HON <- scoresPCATotal_HON %>% 
  # convert long to wide
  tidyr::gather(key, value, -Species) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# loadings
loadingsPCATotal_HON <- as.data.frame(unclass(pcaTotal_HON$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()

loadingsPCATotal_HON <- loadingsPCATotal_HON %>% 
  # convert long to wide
  tidyr::gather(key, value, -trait) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# scalar to adjust arrow size
sc_mi <- 7

loadings_sc_HON <- loadingsPCATotal_HON %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")), funs(.*sc_mi)) %>% 
  # variable names
  mutate(trait = c("Bill width","Hand-Wing Index", "Total culmen", "Tail length", "Tarsus", "Body mass", "Wing length"))

#######################

#### --------------------------------------------------------------
## PCA ## FLORENCIA
#### --------------------------------------------------------------

#table just with columns of morphological traits and species names as rownames
tr_mi_z_FLOR <- FLOR %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)%>% 
  column_to_rownames(var = "Species")

#Run PCA
pcaTotal_FLOR <- princomp(tr_mi_z_FLOR, cor = TRUE, scores = TRUE)
summary(pcaTotal_FLOR)

#Scores
scoresPCATotal_FLOR <- as.data.frame(pcaTotal_FLOR$scores) %>% 
  tibble::rownames_to_column("Species")

scoresPCATotal_FLOR <- scoresPCATotal_FLOR %>% 
  # convert long to wide
  tidyr::gather(key, value, -Species) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# loadings
loadingsPCATotal_FLOR <- as.data.frame(unclass(pcaTotal_FLOR$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()

loadingsPCATotal_FLOR <- loadingsPCATotal_FLOR %>% 
  # convert long to wide
  tidyr::gather(key, value, -trait) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# scalar to adjust arrow size
sc_mi <- 7

loadings_sc_FLOR <- loadingsPCATotal_FLOR %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")), funs(.*sc_mi)) %>% 
  # variable names
  mutate(trait = c("Bill width","Hand-Wing Index", "Total culmen", "Tail length", "Tarsus", "Body mass", "Wing length"))

#######################

#### --------------------------------------------------------------
## PCA ## FUSAGASUGA
#### --------------------------------------------------------------

#table just with columns of morphological traits and species names as rownames
tr_mi_z_FUSA <- FUSA %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)%>% 
  column_to_rownames(var = "Species")

#Run PCA
pcaTotal_FUSA <- princomp(tr_mi_z_FUSA, cor = TRUE, scores = TRUE)
summary(pcaTotal_FUSA)

#Scores
scoresPCATotal_FUSA <- as.data.frame(pcaTotal_FUSA$scores) %>% 
  tibble::rownames_to_column("Species")

scoresPCATotal_FUSA <- scoresPCATotal_FUSA %>% 
  # convert long to wide
  tidyr::gather(key, value, -Species) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# loadings
loadingsPCATotal_FUSA <- as.data.frame(unclass(pcaTotal_FUSA$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()

loadingsPCATotal_FUSA <- loadingsPCATotal_FUSA %>% 
  # convert long to wide
  tidyr::gather(key, value, -trait) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# scalar to adjust arrow size
sc_mi <- 7

loadings_sc_FUSA <- loadingsPCATotal_FUSA %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")), funs(.*sc_mi)) %>% 
  # variable names
  mutate(trait = c("Bill width","Hand-Wing Index", "Total culmen", "Tail length", "Tarsus", "Body mass", "Wing length"))

#######################

#Subset for species in each Region and period

MIR_H <- trait_z %>% 
  filter(MIR_H_corr >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_MIR, by = "Species") 

MIR_M <- trait_z %>% 
  filter(MIR_M >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_MIR, by = "Species")

BARB_H <- trait_z %>% 
  filter(BARB_H_corr >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_BARB, by = "Species") 

BARB_M <- trait_z %>% 
  filter(BARB_M >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_BARB, by = "Species")

SAGU_H <- trait_z %>% 
  filter(SAGU_H_corr >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_SAGU, by = "Species") 

SAGU_M <- trait_z %>% 
  filter(SAGU_M >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_SAGU, by = "Species")

TOCHE_H <- trait_z %>% 
  filter(TOCHE_H_corr >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_TOCHE, by = "Species") 

TOCHE_M <- trait_z %>% 
  filter(TOCHE_M >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_TOCHE, by = "Species")

HON_H <- trait_z %>% 
  filter(HON_H_corr >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_HON, by = "Species") 

HON_M <- trait_z %>% 
  filter(HON_M >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_HON, by = "Species")

FLOR_H <- trait_z %>% 
  filter(FLOR_H_corr >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FLOR, by = "Species") 

FLOR_M <- trait_z %>% 
  filter(FLOR_M >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FLOR, by = "Species")

FUSA_H <- trait_z %>% 
  filter(FUSA_H_corr >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FUSA, by = "Species") 

FUSA_M <- trait_z %>% 
  filter(FUSA_M >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FUSA, by = "Species")

# kernel density estimation for each Region and each period
pc_raw_MIR_H <- MIR_H %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "Species")

pc_raw_MIR_M <- MIR_M %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "Species")

pc_raw_BARB_H <- BARB_H %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "Species")

pc_raw_BARB_M <- BARB_M %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "Species")

pc_raw_SAGU_M <- SAGU_M %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "Species")

pc_raw_SAGU_H <- SAGU_H %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "Species")

pc_raw_TOCHE_H <- TOCHE_H %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "Species")

pc_raw_TOCHE_M <- TOCHE_M %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "Species")

pc_raw_HON_H <- HON_H %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "Species")

pc_raw_HON_M <- HON_M %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "Species")

pc_raw_FLOR_H <- FLOR_H %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "Species")

pc_raw_FLOR_M <- FLOR_M %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "Species")

pc_raw_FUSA_H <- FUSA_H %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "Species")

pc_raw_FUSA_M <- FUSA_M %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "Species")

# optimal bandwidth estimation
hpi_MIR_H <- Hpi(x = pc_raw_MIR_H)
hpi_MIR_M <- Hpi(x = pc_raw_MIR_M)

hpi_BARB_H <- Hpi(x = pc_raw_BARB_H)
hpi_BARB_M <- Hpi(x = pc_raw_BARB_M)

hpi_SAGU_H <- Hpi(x = pc_raw_SAGU_H)
hpi_SAGU_M <- Hpi(x = pc_raw_SAGU_M)

hpi_TOCHE_H <- Hpi(x = pc_raw_TOCHE_H)
hpi_TOCHE_M <- Hpi(x = pc_raw_TOCHE_M)

hpi_HON_H <- Hpi(x = pc_raw_HON_H)
hpi_HON_M <- Hpi(x = pc_raw_HON_M)

hpi_FLOR_H <- Hpi(x = pc_raw_FLOR_H)
hpi_FLOR_M <- Hpi(x = pc_raw_FLOR_M)

hpi_FUSA_H <- Hpi(x = pc_raw_FUSA_H)
hpi_FUSA_M <- Hpi(x = pc_raw_FUSA_M)

# kernel density estimation Miraflores_Historical  
est_MIR_H <- kde(x = pc_raw_MIR_H, H = hpi_MIR_H, compute.cont = TRUE)  

den_MIR_H <- list(est_MIR_H$eval.points[[1]], est_MIR_H$eval.points[[2]], est_MIR_H$estimate)
names(den_MIR_H) <- c("x", "y", "z")
dimnames(den_MIR_H$z) <- list(den_MIR_H$x, den_MIR_H$y)
dcc_MIR_H <- melt(den_MIR_H$z)

# source: kernel function
## --------------------------------------------------------------
## Name: 2.1.3-kernel-function.R
## Description: Function to calculate kernel density probablities
## Date: October 2017
## Outputs: Function named 'cl'
## Args:
## df = dataframe of kernel density data
## prob = Probabilty level e.g. 0.95 (95% confidence level)
## --------------------------------------------------------------

cl <- function(df, prob) {
  dx <- diff(df$x[1:2])
  dy <- diff(df$y[1:2])
  sz <- sort(df$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}

# 0.5 probability kernel
cl_50_MIR_H <- cl(df = den_MIR_H, prob = 0.50)
# 0.95 probability kernel
cl_95_MIR_H <- cl(df = den_MIR_H, prob = 0.95)
# 0.99 probability kernel
cl_99_MIR_H <- cl(df = den_MIR_H, prob = 0.99)

# kernel density estimation miraflores_Modern

est_MIR_M <- kde(x = pc_raw_MIR_M, H = hpi_MIR_M, compute.cont = TRUE)  

den_MIR_M <- list(est_MIR_M$eval.points[[1]], est_MIR_M$eval.points[[2]], est_MIR_M$estimate)
names(den_MIR_M) <- c("x", "y", "z")
dimnames(den_MIR_M$z) <- list(den_MIR_M$x, den_MIR_M$y)
dcc_MIR_M <- melt(den_MIR_M$z)

# 0.5 probability kernel
cl_50_MIR_M <- cl(df = den_MIR_M, prob = 0.50)
# 0.95 probability kernel
cl_95_MIR_M <- cl(df = den_MIR_M, prob = 0.95)
# 0.99 probability kernel
cl_99_MIR_M <- cl(df = den_MIR_M, prob = 0.99)


# kernel density estimation Barbacoas_Historical  
est_BARB_H <- kde(x = pc_raw_BARB_H, H = hpi_BARB_H, compute.cont = TRUE)  

den_BARB_H <- list(est_BARB_H$eval.points[[1]], est_BARB_H$eval.points[[2]], est_BARB_H$estimate)
names(den_BARB_H) <- c("x", "y", "z")
dimnames(den_BARB_H$z) <- list(den_BARB_H$x, den_BARB_H$y)
dcc_BARB_H <- melt(den_BARB_H$z)

# source: kernel function
## --------------------------------------------------------------
## Name: 2.1.3-kernel-function.R
## Description: Function to calculate kernel density probablities
## Date: October 2017
## Outputs: Function named 'cl'
## Args:
## df = dataframe of kernel density data
## prob = Probabilty level e.g. 0.95 (95% confidence level)
## --------------------------------------------------------------

cl <- function(df, prob) {
  dx <- diff(df$x[1:2])
  dy <- diff(df$y[1:2])
  sz <- sort(df$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}

# 0.5 probability kernel
cl_50_BARB_H <- cl(df = den_BARB_H, prob = 0.50)
# 0.95 probability kernel
cl_95_BARB_H <- cl(df = den_BARB_H, prob = 0.95)
# 0.99 probability kernel
cl_99_BARB_H <- cl(df = den_BARB_H, prob = 0.99)

# kernel density estimation Barbacoas_Modern

est_BARB_M <- kde(x = pc_raw_BARB_M, H = hpi_BARB_M, compute.cont = TRUE)  

den_BARB_M <- list(est_BARB_M$eval.points[[1]], est_BARB_M$eval.points[[2]], est_BARB_M$estimate)
names(den_BARB_M) <- c("x", "y", "z")
dimnames(den_BARB_M$z) <- list(den_BARB_M$x, den_BARB_M$y)
dcc_BARB_M <- melt(den_BARB_M$z)

# 0.5 probability kernel
cl_50_BARB_M <- cl(df = den_BARB_M, prob = 0.50)
# 0.95 probability kernel
cl_95_BARB_M <- cl(df = den_BARB_M, prob = 0.95)
# 0.99 probability kernel
cl_99_BARB_M <- cl(df = den_BARB_M, prob = 0.99)

# kernel density estimation San Agustin_Modern

est_SAGU_M <- kde(x = pc_raw_SAGU_M, H = hpi_SAGU_M, compute.cont = TRUE)  

den_SAGU_M <- list(est_SAGU_M$eval.points[[1]], est_SAGU_M$eval.points[[2]], est_SAGU_M$estimate)
names(den_SAGU_M) <- c("x", "y", "z")
dimnames(den_SAGU_M$z) <- list(den_SAGU_M$x, den_SAGU_M$y)
dcc_SAGU_M <- melt(den_SAGU_M$z)

# 0.5 probability kernel
cl_50_SAGU_M <- cl(df = den_SAGU_M, prob = 0.50)
# 0.95 probability kernel
cl_95_SAGU_M <- cl(df = den_SAGU_M, prob = 0.95)
# 0.99 probability kernel
cl_99_SAGU_M <- cl(df = den_SAGU_M, prob = 0.99)

# kernel density estimation San Agustin_Historical

est_SAGU_H <- kde(x = pc_raw_SAGU_H, H = hpi_SAGU_H, compute.cont = TRUE)  

den_SAGU_H <- list(est_SAGU_H$eval.points[[1]], est_SAGU_H$eval.points[[2]], est_SAGU_H$estimate)
names(den_SAGU_H) <- c("x", "y", "z")
dimnames(den_SAGU_H$z) <- list(den_SAGU_H$x, den_SAGU_H$y)
dcc_SAGU_H <- melt(den_SAGU_H$z)

# 0.5 probability kernel
cl_50_SAGU_H <- cl(df = den_SAGU_H, prob = 0.50)
# 0.95 probability kernel
cl_95_SAGU_H <- cl(df = den_SAGU_H, prob = 0.95)
# 0.99 probability kernel
cl_99_SAGU_H <- cl(df = den_SAGU_H, prob = 0.99)

# kernel density estimation Toche_Modern

est_TOCHE_M <- kde(x = pc_raw_TOCHE_M, H = hpi_TOCHE_M, compute.cont = TRUE)  

den_TOCHE_M <- list(est_TOCHE_M$eval.points[[1]], est_TOCHE_M$eval.points[[2]], est_TOCHE_M$estimate)
names(den_TOCHE_M) <- c("x", "y", "z")
dimnames(den_TOCHE_M$z) <- list(den_TOCHE_M$x, den_TOCHE_M$y)
dcc_TOCHE_M <- melt(den_TOCHE_M$z)

# 0.5 probability kernel
cl_50_TOCHE_M <- cl(df = den_TOCHE_M, prob = 0.50)
# 0.95 probability kernel
cl_95_TOCHE_M <- cl(df = den_TOCHE_M, prob = 0.95)
# 0.99 probability kernel
cl_99_TOCHE_M <- cl(df = den_TOCHE_M, prob = 0.99)

# kernel density estimation Toche_HISTORICAL

est_TOCHE_H <- kde(x = pc_raw_TOCHE_H, H = hpi_TOCHE_H, compute.cont = TRUE)  

den_TOCHE_H <- list(est_TOCHE_H$eval.points[[1]], est_TOCHE_H$eval.points[[2]], est_TOCHE_H$estimate)
names(den_TOCHE_H) <- c("x", "y", "z")
dimnames(den_TOCHE_H$z) <- list(den_TOCHE_H$x, den_TOCHE_H$y)
dcc_TOCHE_H <- melt(den_TOCHE_H$z)

# 0.5 probability kernel
cl_50_TOCHE_H <- cl(df = den_TOCHE_H, prob = 0.50)
# 0.95 probability kernel
cl_95_TOCHE_H <- cl(df = den_TOCHE_H, prob = 0.95)
# 0.99 probability kernel
cl_99_TOCHE_H <- cl(df = den_TOCHE_H, prob = 0.99)

# kernel density estimation Honda_Modern

est_HON_M <- kde(x = pc_raw_HON_M, H = hpi_HON_M, compute.cont = TRUE)  

den_HON_M <- list(est_HON_M$eval.points[[1]], est_HON_M$eval.points[[2]], est_HON_M$estimate)
names(den_HON_M) <- c("x", "y", "z")
dimnames(den_HON_M$z) <- list(den_HON_M$x, den_HON_M$y)
dcc_HON_M <- melt(den_HON_M$z)

# 0.5 probability kernel
cl_50_HON_M <- cl(df = den_HON_M, prob = 0.50)
# 0.95 probability kernel
cl_95_HON_M <- cl(df = den_HON_M, prob = 0.95)
# 0.99 probability kernel
cl_99_HON_M <- cl(df = den_HON_M, prob = 0.99)

# kernel density estimation Honda_Historical

est_HON_H <- kde(x = pc_raw_HON_H, H = hpi_HON_H, compute.cont = TRUE)  

den_HON_H <- list(est_HON_H$eval.points[[1]], est_HON_H$eval.points[[2]], est_HON_H$estimate)
names(den_HON_H) <- c("x", "y", "z")
dimnames(den_HON_H$z) <- list(den_HON_H$x, den_HON_H$y)
dcc_HON_H <- melt(den_HON_H$z)

# 0.5 probability kernel
cl_50_HON_H <- cl(df = den_HON_H, prob = 0.50)
# 0.95 probability kernel
cl_95_HON_H <- cl(df = den_HON_H, prob = 0.95)
# 0.99 probability kernel
cl_99_HON_H <- cl(df = den_HON_H, prob = 0.99)

# kernel density estimation Florencia_Modern

est_FLOR_M <- kde(x = pc_raw_FLOR_M, H = hpi_FLOR_M, compute.cont = TRUE)  

den_FLOR_M <- list(est_FLOR_M$eval.points[[1]], est_FLOR_M$eval.points[[2]], est_FLOR_M$estimate)
names(den_FLOR_M) <- c("x", "y", "z")
dimnames(den_FLOR_M$z) <- list(den_FLOR_M$x, den_FLOR_M$y)
dcc_FLOR_M <- melt(den_FLOR_M$z)

# 0.5 probability kernel
cl_50_FLOR_M <- cl(df = den_FLOR_M, prob = 0.50)
# 0.95 probability kernel
cl_95_FLOR_M <- cl(df = den_FLOR_M, prob = 0.95)
# 0.99 probability kernel
cl_99_FLOR_M <- cl(df = den_FLOR_M, prob = 0.99)

# kernel density estimation Florencia_Historical

est_FLOR_H <- kde(x = pc_raw_FLOR_H, H = hpi_FLOR_H, compute.cont = TRUE)  

den_FLOR_H <- list(est_FLOR_H$eval.points[[1]], est_FLOR_H$eval.points[[2]], est_FLOR_H$estimate)
names(den_FLOR_H) <- c("x", "y", "z")
dimnames(den_FLOR_H$z) <- list(den_FLOR_H$x, den_FLOR_H$y)
dcc_FLOR_H <- melt(den_FLOR_H$z)

# 0.5 probability kernel
cl_50_FLOR_H <- cl(df = den_FLOR_H, prob = 0.50)
# 0.95 probability kernel
cl_95_FLOR_H <- cl(df = den_FLOR_H, prob = 0.95)
# 0.99 probability kernel
cl_99_FLOR_H <- cl(df = den_FLOR_H, prob = 0.99)

# kernel density estimation FUSAGASUGA_Modern

est_FUSA_M <- kde(x = pc_raw_FUSA_M, H = hpi_FUSA_M, compute.cont = TRUE)  

den_FUSA_M <- list(est_FUSA_M$eval.points[[1]], est_FUSA_M$eval.points[[2]], est_FUSA_M$estimate)
names(den_FUSA_M) <- c("x", "y", "z")
dimnames(den_FUSA_M$z) <- list(den_FUSA_M$x, den_FUSA_M$y)
dcc_FUSA_M <- melt(den_FUSA_M$z)

# 0.5 probability kernel
cl_50_FUSA_M <- cl(df = den_FUSA_M, prob = 0.50)
# 0.95 probability kernel
cl_95_FUSA_M <- cl(df = den_FUSA_M, prob = 0.95)
# 0.99 probability kernel
cl_99_FUSA_M <- cl(df = den_FUSA_M, prob = 0.99)

# kernel density estimation FUSAGASUGA_Historical

est_FUSA_H <- kde(x = pc_raw_FUSA_H, H = hpi_FUSA_H, compute.cont = TRUE)  

den_FUSA_H <- list(est_FUSA_H$eval.points[[1]], est_FUSA_H$eval.points[[2]], est_FUSA_H$estimate)
names(den_FUSA_H) <- c("x", "y", "z")
dimnames(den_FUSA_H$z) <- list(den_FUSA_H$x, den_FUSA_H$y)
dcc_FUSA_H <- melt(den_FUSA_H$z)

# 0.5 probability kernel
cl_50_FUSA_H <- cl(df = den_FUSA_H, prob = 0.50)
# 0.95 probability kernel
cl_95_FUSA_H <- cl(df = den_FUSA_H, prob = 0.95)
# 0.99 probability kernel
cl_99_FUSA_H <- cl(df = den_FUSA_H, prob = 0.99)

# save principal component data

PCA_MIR <- trait_z %>% 
  filter(MIR_H_corr >= 0) %>% 
  select(Species,MIR_H_corr, MIR_M,MIR_EXT,MIR_NEW, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_MIR, by = "Species") 

PCA_BARB <- trait_z %>% 
  filter(BARB_H_corr >= 0) %>% 
  select(Species,BARB_H_corr, BARB_M,BARB_EXT,BARB_NEW, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>%  
  left_join(scoresPCATotal_BARB, by = "Species") 

PCA_SAGU <- trait_z %>% 
  filter(SAGU_H_corr >= 0) %>% 
  select(Species,SAGU_H_corr, SAGU_M,SAGU_EXT,SAGU_NEW, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_SAGU, by = "Species") 

PCA_TOCHE <- trait_z %>% 
  filter(TOCHE_H_corr >= 0) %>% 
  select(Species,TOCHE_H_corr, TOCHE_M,TOCHE_EXT,TOCHE_NEW, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_TOCHE, by = "Species") 

PCA_HON <- trait_z %>% 
  filter(HON_H_corr >= 0) %>% 
  select(Species,HON_H_corr,HON_M,HON_EXT, HON_NEW, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_HON, by = "Species") 

PCA_FLOR <- trait_z %>% 
  filter(FLOR_H_corr >= 0) %>% 
  select(Species,FLOR_H_corr, FLOR_M,FLOR_EXT,FLOR_NEW, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FLOR, by = "Species") 

PCA_FUSA <- trait_z %>% 
  filter(FUSA_H_corr >= 0) %>% 
  select(Species,FUSA_H_corr,FUSA_M,FUSA_EXT,FUSA_NEW, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FUSA, by = "Species") 

write.csv(PCA_MIR, file = "PCA_MIR_corr_FINAL.csv", row.names = FALSE)
write.csv(PCA_BARB, file = "PCA_BARB_corr_FINAL.csv", row.names = FALSE)
write.csv(PCA_SAGU, file = "PCA_SAGU_corr_FINAL.csv", row.names = FALSE)
write.csv(PCA_TOCHE, file = "PCA_TOCHE_corr_FINAL.csv", row.names = FALSE)
write.csv(PCA_HON, file = "PCA_HON_corr_FINAL.csv", row.names = FALSE)
write.csv(PCA_FLOR, file = "PCA_FLOR_corr_FINAL.csv", row.names = FALSE)
write.csv(PCA_FUSA, file = "PCA_FUSA_corr_FINAL.csv", row.names = FALSE)

write.csv(scoresPCATotal_MIR, file = "PCA_Total_MIR_corr_FINAL.csv")
write.csv(scoresPCATotal_BARB, file = "PCA_Total_BARB_corr_FINAL.csv")
write.csv(scoresPCATotal_TOCHE, file = "PCA_Total_TOCHE_corr_FINAL.csv")
write.csv(scoresPCATotal_SAGU, file = "PCA_Total_SAGU_corr_FINAL.csv")
write.csv(scoresPCATotal_HON, file = "PCA_Total_HON_corr_FINAL.csv")
write.csv(scoresPCATotal_FLOR, file = "PCA_Total_FLOR_corr_FINAL.csv")
write.csv(scoresPCATotal_FUSA, file = "PCA_Total_FUSA_corr_FINAL.csv")

write.csv(loadingsPCATotal_MIR, file = "Loadings_Total_MIR_corr_FINAL.csv", row.names = FALSE)
write.csv(loadingsPCATotal_BARB, file = "Loadings_Total_BARB_corr_FINAL.csv", row.names = FALSE)
write.csv(loadingsPCATotal_SAGU, file = "Loadings_Total_SAGU.csv_corr_FINAL", row.names = FALSE)
write.csv(loadingsPCATotal_TOCHE, file = "Loadings_Total_TOCHE_corr_FINAL.csv", row.names = FALSE)
write.csv(loadingsPCATotal_HON, file = "Loadings_Total_HON_corr_FINAL.csv", row.names = FALSE)
write.csv(loadingsPCATotal_FLOR, file = "Loadings_Total_FLOR_corr_FINAL.csv", row.names = FALSE)
write.csv(loadingsPCATotal_FUSA, file = "Loadings_Total_FUSA_corr_FINAL.csv", row.names = FALSE)

#### PCA plots
library(ggplot2)

# colour palette
col_pal <- colorRampPalette(c("grey30","grey60", "grey80", "white"))(200)


### Folowing lines indentify species that correspond to extinction,
#recolonization or new additions in different periods

extir_MIR <- subset(PCA_MIR, MIR_EXT == 1)
ne_MIR <- subset(PCA_MIR, MIR_NEW == 1)

extir_BARB <- subset(PCA_BARB, BARB_EXT == 1)
ne_BARB <- subset(PCA_BARB, BARB_NEW == 1)

extir_SAGU <- subset(PCA_SAGU, SAGU_EXT == 1)
ne_SAGU <- subset(PCA_SAGU, SAGU_NEW == 1)

extir_TOCHE <- subset(PCA_TOCHE, TOCHE_EXT == 1)
ne_TOCHE <- subset(PCA_TOCHE, TOCHE_NEW == 1)

extir_HON <- subset(PCA_HON, HON_EXT == 1)
ne_HON <- subset(PCA_HON, HON_NEW == 1)

extir_FLOR <- subset(PCA_FLOR, FLOR_EXT == 1)
ne_FLOR <- subset(PCA_FLOR, FLOR_NEW == 1)

extir_FUSA <- subset(PCA_FUSA, FUSA_EXT == 1)
ne_FUSA <- subset(PCA_FUSA, FUSA_NEW == 1)

# plot Miraflores Historical
pca_plot_MIR_H <- ggplot(dcc_MIR_H, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_MIR, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_MIR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = ne_MIR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "light blue") +
  
  # add arrows
  geom_segment(data = loadings_sc_MIR, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_MIR_H, colour = "grey30", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_95_MIR_H, colour = "grey60", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_99_MIR_H, colour = "grey70", linewidth = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_MIR, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_MIR, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_MIR, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (67%)", y = "PC2 - Dispersal ability (20%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Miraflores Historical", color="black", size = 5)

# display plot
windows()
pca_plot_MIR_H

# plot MIRAFLORES Modern
pca_plot_MIR_M <- ggplot(dcc_MIR_M, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_MIR, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_MIR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = ne_MIR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "light blue") +
  
  # add arrows
  geom_segment(data = loadings_sc_MIR, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_MIR_M, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_MIR_M, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_MIR_M, colour = "grey70", size = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_MIR, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_MIR, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_MIR, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (67%)", y = "PC2 - Dispersal ability (20%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Miraflores Modern", color="black", size = 5)

# display plot
windows()
pca_plot_MIR_M


# plot Barbacoas Historical
pca_plot_BARB_H <- ggplot(dcc_BARB_H, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_BARB, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_BARB, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = ne_BARB, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "light blue") +
  
  # add arrows
  geom_segment(data = loadings_sc_BARB, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_BARB_H, colour = "grey30", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_95_BARB_H, colour = "grey60", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_99_BARB_H, colour = "grey70", linewidth = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_BARB, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_BARB, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_BARB, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (63%)", y = "PC2 - Dispersal ability (19%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Barbacoas Historical", color="black", size = 5)

# display plot
windows()
pca_plot_BARB_H

# plot Barbacoas Modern
pca_plot_BARB_M <- ggplot(dcc_BARB_M, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_BARB, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_BARB, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = ne_BARB, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "light blue") +
  
  # add arrows
  geom_segment(data = loadings_sc_BARB, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_BARB_M, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_BARB_M, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_BARB_M, colour = "grey70", size = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_BARB, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_BARB, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_BARB, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (63%)", y = "PC2 - Dispersal ability (19%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Barbacoas Modern", color="black", size = 5)

# display plot
windows()
pca_plot_BARB_M

# plot San Agustin Historical
pca_plot_SAGU_H <- ggplot(dcc_SAGU_H, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_SAGU, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_SAGU, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = ne_SAGU, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "light blue") +
  
  # add arrows
  geom_segment(data = loadings_sc_SAGU, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_SAGU_H, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_SAGU_H, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_SAGU_H, colour = "grey70", size = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_SAGU, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_SAGU, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_SAGU, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (65%)", y = "PC2 - Dispersal ability (20%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="San Agustin Historical", color="black", size = 5)

# display plot
windows()
pca_plot_SAGU_H

# plot San Agustin Modern
pca_plot_SAGU_M <- ggplot(dcc_SAGU_M, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_SAGU, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_SAGU, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = ne_SAGU, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "light blue") +
  
  # add arrows
  geom_segment(data = loadings_sc_SAGU, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_SAGU_M, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_SAGU_M, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_SAGU_M, colour = "grey70", size = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_SAGU, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_SAGU, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_SAGU, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (65%)", y = "PC2 - Dispersal ability (20%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="San Agustin Modern", color="black", size = 5)

# display plot
windows()
pca_plot_SAGU_M

# plot Toche Historical
pca_plot_TOCHE_H <- ggplot(dcc_TOCHE_H, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_TOCHE, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_TOCHE, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = ne_TOCHE, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "light blue") +
  
  # add arrows
  geom_segment(data = loadings_sc_TOCHE, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_TOCHE_H, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_TOCHE_H, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_TOCHE_H, colour = "grey70", size = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_TOCHE, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_TOCHE, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_TOCHE, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (66%)", y = "PC2 - Dispersal ability (20%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Toche Historical", color="black", size = 5)

# display plot
windows()
pca_plot_TOCHE_H

# plot Toche Modern
pca_plot_TOCHE_M <- ggplot(dcc_TOCHE_M, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_TOCHE, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_TOCHE, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = ne_TOCHE, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "light blue") +
  
  # add arrows
  geom_segment(data = loadings_sc_TOCHE, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_TOCHE_M, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_TOCHE_M, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_TOCHE_M, colour = "grey70", size = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_TOCHE, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_TOCHE, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_TOCHE, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (65%)", y = "PC2 - Dispersal ability (20%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Toche Modern", color="black", size = 5)

# display plot
windows()
pca_plot_TOCHE_M


# plot Honda Historical
pca_plot_HON_H <- ggplot(dcc_HON_H, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_HON, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_HON, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = ne_HON, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "light blue") +
  
  # add arrows
  geom_segment(data = loadings_sc_HON, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_HON_H, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_HON_H, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_HON_H, colour = "grey70", size = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_HON, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_HON, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_HON, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (64%)", y = "PC2 - Dispersal ability (19%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Honda Historical", color="black", size = 5)

# display plot
windows()
pca_plot_HON_H

# plot Honda Modern
pca_plot_HON_M <- ggplot(dcc_HON_M, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_HON, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_HON, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = ne_HON, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "light blue") +
  
  # add arrows
  geom_segment(data = loadings_sc_HON, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_HON_M, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_HON_M, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_HON_M, colour = "grey70", size = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_HON, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_HON, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_HON, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (64%)", y = "PC2 - Dispersal ability (19%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Honda Modern", color="black", size = 5)

# display plot
windows()
pca_plot_HON_M


# plot Florencia Historical
pca_plot_FLOR_H <- ggplot(dcc_FLOR_H, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_FLOR, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_FLOR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = ne_FLOR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "light blue") +
  
  # add arrows
  geom_segment(data = loadings_sc_FLOR, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_FLOR_H, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_FLOR_H, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_FLOR_H, colour = "grey70", size = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_FLOR, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_FLOR, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_FLOR, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (65%)", y = "PC2 - Dispersal ability (18%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Florencia Historical", color="black", size = 5)

# display plot
windows()
pca_plot_FLOR_H

# plot Florencia Modern
pca_plot_FLOR_M <- ggplot(dcc_FLOR_M, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_FLOR, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_FLOR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = ne_FLOR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "light blue") +
  
  # add arrows
  geom_segment(data = loadings_sc_FLOR, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_FLOR_M, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_FLOR_M, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_FLOR_M, colour = "grey70", size = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_FLOR, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_FLOR, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_FLOR, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (65%)", y = "PC2 - Dispersal ability (18%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Florencia Modern", color="black", size = 5)

# display plot
windows()
pca_plot_FLOR_M

# plot Fusagasuga Historical
pca_plot_FUSA_H <- ggplot(dcc_FUSA_H, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_FUSA, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_FUSA, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = ne_FUSA, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "light blue") +
  
  # add arrows
  geom_segment(data = loadings_sc_FUSA, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_FUSA_H, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_FUSA_H, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_FUSA_H, colour = "grey70", size = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_FUSA, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_FUSA, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_FUSA, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (63%)", y = "PC2 - Dispersal ability (21%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Fusagasuga Historical", color="black", size = 5)

# display plot
windows()
pca_plot_FUSA_H

# plot Fusagasuga Modern
pca_plot_FUSA_M <- ggplot(dcc_FUSA_M, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_FUSA, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_FUSA, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = ne_FUSA, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "light blue") +
  
  # add arrows
  geom_segment(data = loadings_sc_FUSA, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_FUSA_M, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_FUSA_M, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_FUSA_M, colour = "grey70", size = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_FUSA, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_FUSA, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_FUSA, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (63%)", y = "PC2 - Dispersal ability (22%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Fusagasuga Modern", color="black", size = 5)

# display plot
windows()
pca_plot_FUSA_M

#################HYPER VOLUMES##################################


#### Estimating hypervolumes by period for all regions (hypervolume overlap) using PCA

# dataframes of trait data (PCA scores)

#MIRAFLORES
PCA_MIR_H_A <- subset(PCA_MIR, MIR_H_corr >= 1)

PCA_MIR_H_A <-  PCA_MIR_H_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3, Comp.4, Comp.5) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

PCA_MIR_M_A <- subset(PCA_MIR, MIR_M >= 1)

PCA_MIR_M_A <-  PCA_MIR_M_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3, Comp.4, Comp.5) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

#BARBACOAS
PCA_BARB_H_A <- subset(PCA_BARB, BARB_H_corr >= 1)
  
PCA_BARB_H_A <-  PCA_BARB_H_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3, Comp.4, Comp.5) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

PCA_BARB_M_A <- subset(PCA_BARB, BARB_M >= 1)

PCA_BARB_M_A <-  PCA_BARB_M_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3, Comp.4, Comp.5) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

#TOCHE
PCA_TOCHE_H_A <- subset(PCA_TOCHE, TOCHE_H_corr >= 1)

PCA_TOCHE_H_A <-  PCA_TOCHE_H_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3, Comp.4, Comp.5) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species") 

PCA_TOCHE_M_A <- subset(PCA_TOCHE, TOCHE_M >= 1)

PCA_TOCHE_M_A <-  PCA_TOCHE_M_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3, Comp.4, Comp.5) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

#SAN AGUSTIN
PCA_SAGU_H_A <- subset(PCA_SAGU, SAGU_H_corr >= 1)

PCA_SAGU_H_A <-  PCA_SAGU_H_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3, Comp.4, Comp.5) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species") 

PCA_SAGU_M_A <- subset(PCA_SAGU, SAGU_M >= 1)

PCA_SAGU_M_A <-  PCA_SAGU_M_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3, Comp.4, Comp.5) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

#HONDA
PCA_HON_H_A <- subset(PCA_HON, HON_H_corr >= 1)

PCA_HON_H_A <-  PCA_HON_H_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3, Comp.4, Comp.5) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species") 

PCA_HON_M_A <- subset(PCA_HON, HON_M >= 1)

PCA_HON_M_A <-  PCA_HON_M_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3, Comp.4, Comp.5) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

#FLORENCIA
PCA_FLOR_H_A <- subset(PCA_FLOR, FLOR_H_corr >= 1)

PCA_FLOR_H_A <-  PCA_FLOR_H_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3, Comp.4, Comp.5) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species") 

PCA_FLOR_M_A <- subset(PCA_FLOR, FLOR_M >= 1)

PCA_FLOR_M_A <-  PCA_FLOR_M_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3, Comp.4, Comp.5) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

#FUSAGASUGA
PCA_FUSA_H_A <- subset(PCA_FUSA, FUSA_H_corr >= 1)

PCA_FUSA_H_A <-  PCA_FUSA_H_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3, Comp.4, Comp.5)  %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

PCA_FUSA_M_A <- subset(PCA_FUSA, FUSA_M >= 1)

PCA_FUSA_M_A <-  PCA_FUSA_M_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3, Comp.4, Comp.5) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")


# hypervolumes: Takes a long time to run if all components included. To test, first run using just 3 components by changing code above to:
# PCA_BARB_H_A <-  PCA_BARB_H_A %>% 
# select(Species, Comp.1, Comp.2, Comp.3) %>% 
# remove_rownames %>% 
# column_to_rownames(var = "Species")


library(hypervolume)

#Miraflores
set.seed(3)
MIR_H_hyper_mi <- hypervolume::hypervolume_svm(PCA_MIR_H_A, name = "MIR_H")

saveRDS(MIR_H_hyper_mi, "MIR_H_hyper_mi5_corr_FINAL.rds")

set.seed(3)
MIR_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_MIR_M_A, name = "MIR_M")

saveRDS(MIR_M_hyper_mi, "MIR_M_hyper_mi5_corr_FINAL.rds")

# load: hyper volumes
#MIR_H_hyper_mi <- readRDS("MIR_H_hyper_mi5.rds")
#MIR_M_hyper_mi <- readRDS("MIR_M_hyper_mi5.rds")


# set hypervolume comparisons 

set.seed(3)
bm_set_mi_MIRHvsMIRM <- hypervolume::hypervolume_set(MIR_H_hyper_mi, MIR_M_hyper_mi, check.memory = FALSE)

# overlap statistics
bm_over_mi <- hypervolume::hypervolume_overlap_statistics(bm_set_mi_MIRHvsMIRM)

# summarise volumes Historical  vs Modern

hypervolume::get_volume(bm_set_mi_MIRHvsMIRM)

# Plot 3D hypervolumes 
library(rgl)
library(alphahull)

#Historical vs Modern
bm_set_mi_MIRHvsMIRM@HVList$HV1 <- NULL
bm_set_mi_MIRHvsMIRM@HVList$HV2 <- NULL
bm_set_mi_MIRHvsMIRM@HVList$Union <- NULL

# trait names
names <- c("PC1", "PC2","PC3", "PC4", "PC5")

colnames(bm_set_mi_MIRHvsMIRM@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_MIRHvsMIRM@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_MIRHvsMIRM@HVList$Unique_2@RandomPoints) <- names

# plot hypervolumes 
hypervolume::plot.HypervolumeList(bm_set_mi_MIRHvsMIRM, show.3d = TRUE, plot.3d.axes.id = c(1,2,3), show.random = FALSE, show.data = TRUE, show.centroid = FALSE, show.density = TRUE, cex.random = 5, colors = c("azure3", "red","cornflowerblue"))

#Fancy plot 3D Hypervolume
######################################
#Data frame for Plotly Hyper Volume Plot
inter <- data.frame(bm_set_mi_MIRHvsMIRM@HVList$Intersection@RandomPoints)
inter <- inter %>% 
  mutate(section = "intersection")

uni1 <- data.frame(bm_set_mi_MIRHvsMIRM@HVList$Unique_1@RandomPoints) 
uni1 <- uni1 %>% 
  mutate(section = "Unique Historical")

uni2 <- data.frame(bm_set_mi_MIRHvsMIRM@HVList$Unique_2@RandomPoints)
uni2 <- uni2 %>% 
  mutate(section = "Unique Modern")

HV <- rbind(inter, uni1, uni2)
write.csv(HV,"HyperVolume_MIR.csv")

colnames(bm_set_mi_MIRHvsMIRM@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_MIRHvsMIRM@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_MIRHvsMIRM@HVList$Unique_2@RandomPoints) <- names


library(plotly)

HV <- read.csv("HyperVolume_MIR.csv", h = T)


pal <- c("azure3","red", "cornflowerblue" )


fig <- plot_ly(HV, x = ~PC1, y = ~(-PC2), z = ~(-PC3), type = 'scatter3d', mode = 'markers', color = ~section, colors = pal, marker = list(size = 3,opacity = 0.6, line = list(color = 'grey66', width = 2)))

fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                   yaxis = list(title = 'PC2'),
                                   zaxis = list(title = 'PC3')))

fig

######################################

#Barbacoas 
set.seed(3)
BARB_H_hyper_mi <- hypervolume::hypervolume_svm(PCA_BARB_H_A, name = "BARB_H")

saveRDS(BARB_H_hyper_mi, "BARB_H_hyper_mi5_corr_FINAL.rds")

set.seed(3)
BARB_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_BARB_M_A, name = "BARB_M")

saveRDS(BARB_M_hyper_mi, "BARB_M_hyper_mi5_corr_FINAL.rds")

# load: hyper volumes
#BARB_H_hyper_mi <- readRDS("BARB_H_hyper_mi5.rds")
#BARB_M_hyper_mi <- readRDS("BARB_M_hyper_mi5.rds")


# set hypervolume comparisons 

set.seed(3)
bm_set_mi_BARBHvsBARBM <- hypervolume::hypervolume_set(BARB_H_hyper_mi, BARB_M_hyper_mi, check.memory = FALSE)

# overlap statistics
bm_over_mi <- hypervolume::hypervolume_overlap_statistics(bm_set_mi_BARBHvsBARBM)

# summarise volumes Historical  vs Modern

hypervolume::get_volume(bm_set_mi_BARBHvsBARBM)

# Plot 3D hypervolumes 
library(rgl)
library(alphahull)

#Historical vs Modern
bm_set_mi_BARBHvsBARBM@HVList$HV1 <- NULL
bm_set_mi_BARBHvsBARBM@HVList$HV2 <- NULL
bm_set_mi_BARBHvsBARBM@HVList$Union <- NULL

# trait names
names <- c("PC1", "PC2","PC3", "PC4", "PC5")

colnames(bm_set_mi_BARBHvsBARBM@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_BARBHvsBARBM@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_BARBHvsBARBM@HVList$Unique_2@RandomPoints) <- names

# plot hypervolumes 
hypervolume::plot.HypervolumeList(bm_set_mi_BARBHvsBARBM, show.3d = TRUE, plot.3d.axes.id = c(1,2,3), show.random = FALSE, show.data = TRUE, show.centroid = FALSE, show.density = TRUE, cex.random = 5, colors = c("azure3", "red","cornflowerblue"))

#Fancy plot 3D Hypervolume
######################################
#Data frame for Plotly Hyper Volume Plot
inter <- data.frame(bm_set_mi_BARBHvsBARBM@HVList$Intersection@RandomPoints)
inter <- inter %>% 
  mutate(section = "intersection")

uni1 <- data.frame(bm_set_mi_BARBHvsBARBM@HVList$Unique_1@RandomPoints) 
uni1 <- uni1 %>% 
  mutate(section = "Unique Historical")

uni2 <- data.frame(bm_set_mi_BARBHvsBARBM@HVList$Unique_2@RandomPoints)
uni2 <- uni2 %>% 
  mutate(section = "Unique Modern")

HV <- rbind(inter, uni1, uni2)
write.csv(HV,"HyperVolume_BARB.csv")

colnames(bm_set_mi_BARBHvsBARBM@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_BARBHvsBARBM@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_BARBHvsBARBM@HVList$Unique_2@RandomPoints) <- names


library(plotly)

HV <- read.csv("HyperVolume_BARB.csv", h = T)


pal <- c("azure3","red", "cornflowerblue" )


fig <- plot_ly(HV, x = ~PC1, y = ~(-PC2), z = ~(-PC3), type = 'scatter3d', mode = 'markers', color = ~section, colors = pal, marker = list(size = 3,opacity = 0.6, line = list(color = 'grey66', width = 2)))

fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                   yaxis = list(title = 'PC2'),
                                   zaxis = list(title = 'PC3')))

fig

######################################

#TOCHE
set.seed(3)
TOCHE_H_hyper_mi <- hypervolume::hypervolume_svm(PCA_TOCHE_H_A, name = "TOCHE_H")

saveRDS(TOCHE_H_hyper_mi, "TOCHE_H_hyper_mi5_corr_FINAL.rds")

set.seed(3)
TOCHE_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_TOCHE_M_A, name = "TOCHE_M")

saveRDS(TOCHE_M_hyper_mi, "TOCHE_M_hyper_mi5_corr_FINAL.rds")

# load: hyper volumes
#TOCHE_H_hyper_mi <- readRDS("TOCHE_H_hyper_mi5.rds")
#TOCHE_M_hyper_mi <- readRDS("TOCHE_M_hyper_mi5.rds")


# set hypervolume comparisons 

set.seed(3)
bm_set_mi_TOCHEHvsTOCHEM <- hypervolume::hypervolume_set(TOCHE_H_hyper_mi, TOCHE_M_hyper_mi, check.memory = FALSE)

# overlap statistics
bm_over_mi <- hypervolume::hypervolume_overlap_statistics(bm_set_mi_TOCHEHvsTOCHEM)

# summarise volumes Historical  vs Modern

hypervolume::get_volume(bm_set_mi_TOCHEHvsTOCHEM)

# Plot 3D hypervolumes 
#Historical vs Modern
bm_set_mi_TOCHEHvsTOCHEM@HVList$HV1 <- NULL
bm_set_mi_TOCHEHvsTOCHEM@HVList$HV2 <- NULL
bm_set_mi_TOCHEHvsTOCHEM@HVList$Union <- NULL

# trait names
names <- c("PC1", "PC2","PC3", "PC4", "PC5")

colnames(bm_set_mi_TOCHEHvsTOCHEM@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_TOCHEHvsTOCHEM@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_TOCHEHvsTOCHEM@HVList$Unique_2@RandomPoints) <- names

# plot hypervolumes 
hypervolume::plot.HypervolumeList(bm_set_mi_TOCHEHvsTOCHEM, show.3d = TRUE, plot.3d.axes.id = c(1,2,3), show.random = FALSE, show.data = TRUE, show.centroid = FALSE, show.density = TRUE, cex.random = 5, colors = c("azure3", "red","cornflowerblue"))

#Fancy plot 3D Hypervolume
######################################
#Data frame for Plotly Hyper Volume Plot
inter <- data.frame(bm_set_mi_TOCHEHvsTOCHEM@HVList$Intersection@RandomPoints)
inter <- inter %>% 
  mutate(section = "intersection")

uni1 <- data.frame(bm_set_mi_TOCHEHvsTOCHEM@HVList$Unique_1@RandomPoints) 
uni1 <- uni1 %>% 
  mutate(section = "Unique Historical")

uni2 <- data.frame(bm_set_mi_TOCHEHvsTOCHEM@HVList$Unique_2@RandomPoints)
uni2 <- uni2 %>% 
  mutate(section = "Unique Modern")

HV <- rbind(inter, uni1, uni2)
write.csv(HV,"HyperVolume_TOCHE.csv")

colnames(bm_set_mi_TOCHEHvsTOCHEM@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_TOCHEHvsTOCHEM@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_TOCHEHvsTOCHEM@HVList$Unique_2@RandomPoints) <- names


library(plotly)

HV <- read.csv("HyperVolume_TOCHE.csv", h = T)

pal <- c("azure3","red", "cornflowerblue" )


fig <- plot_ly(HV, x = ~PC1, y = ~(-PC2), z = ~(-PC3), type = 'scatter3d', mode = 'markers', color = ~section, colors = pal, marker = list(size = 3,opacity = 0.6, line = list(color = 'grey66', width = 2)))

fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                   yaxis = list(title = 'PC2'),
                                   zaxis = list(title = 'PC3')))

fig

#SAN AGUSTIN
set.seed(3)
SAGU_H_hyper_mi <- hypervolume::hypervolume_svm(PCA_SAGU_H_A, name = "SAGU_H")

saveRDS(SAGU_H_hyper_mi, "SAGU_H_hyper_mi5_corr_FINAL.rds")

set.seed(3)
SAGU_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_SAGU_M_A, name = "SAGU_M")

saveRDS(SAGU_M_hyper_mi, "SAGU_M_hyper_mi5_corr_FINAL.rds")

# load: hyper volumes
#SAGU_H_hyper_mi <- readRDS("SAGU_H_hyper_mi5.rds")
#SAGU_M_hyper_mi <- readRDS("SAGU_M_hyper_mi5.rds")


# set hypervolume comparisons 

set.seed(3)
bm_set_mi_SAGUHvsSAGUM <- hypervolume::hypervolume_set(SAGU_H_hyper_mi, SAGU_M_hyper_mi, check.memory = FALSE)

# overlap statistics
bm_over_mi <- hypervolume::hypervolume_overlap_statistics(bm_set_mi_SAGUHvsSAGUM)

# summarise volumes Historical  vs Modern

hypervolume::get_volume(bm_set_mi_SAGUHvsSAGUM)

# Plot 3D hypervolumes 
#Historical vs Modern
bm_set_mi_SAGUHvsSAGUM@HVList$HV1 <- NULL
bm_set_mi_SAGUHvsSAGUM@HVList$HV2 <- NULL
bm_set_mi_SAGUHvsSAGUM@HVList$Union <- NULL

# trait names
names <- c("PC1", "PC2","PC3", "PC4", "PC5")

colnames(bm_set_mi_SAGUHvsSAGUM@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_SAGUHvsSAGUM@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_SAGUHvsSAGUM@HVList$Unique_2@RandomPoints) <- names

# plot hypervolumes 
hypervolume::plot.HypervolumeList(bm_set_mi_SAGUHvsSAGUM, show.3d = TRUE, plot.3d.axes.id = c(1,2,3), show.random = FALSE, show.data = TRUE, show.centroid = FALSE, show.density = TRUE, cex.random = 5, colors = c("azure3", "red","cornflowerblue"))

#Fancy plot 3D Hypervolume
######################################
#Data frame for Plotly Hyper Volume Plot
inter <- data.frame(bm_set_mi_SAGUHvsSAGUM@HVList$Intersection@RandomPoints)
inter <- inter %>% 
  mutate(section = "intersection")

uni1 <- data.frame(bm_set_mi_SAGUHvsSAGUM@HVList$Unique_1@RandomPoints) 
uni1 <- uni1 %>% 
  mutate(section = "Unique Historical")

uni2 <- data.frame(bm_set_mi_SAGUHvsSAGUM@HVList$Unique_2@RandomPoints)
uni2 <- uni2 %>% 
  mutate(section = "Unique Modern")

HV <- rbind(inter, uni1, uni2)
write.csv(HV,"HyperVolume_SAGU.csv")

colnames(bm_set_mi_SAGUHvsSAGUM@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_SAGUHvsSAGUM@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_SAGUHvsSAGUM@HVList$Unique_2@RandomPoints) <- names


library(plotly)

HV <- read.csv("HyperVolume_SAGU.csv", h = T)

pal <- c("azure3","red", "cornflowerblue" )


fig <- plot_ly(HV, x = ~PC1, y = ~(-PC2), z = ~(-PC3), type = 'scatter3d', mode = 'markers', color = ~section, colors = pal, marker = list(size = 3,opacity = 0.6, line = list(color = 'grey66', width = 2)))

fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                   yaxis = list(title = 'PC2'),
                                   zaxis = list(title = 'PC3')))

fig


#HONDA
set.seed(3)
HON_H_hyper_mi <- hypervolume::hypervolume_svm(PCA_HON_H_A, name = "HON_H")

saveRDS(HON_H_hyper_mi, "HON_H_hyper_mi5_corr_FINAL.rds")

set.seed(3)
HON_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_HON_M_A, name = "HON_M")

saveRDS(HON_M_hyper_mi, "HON_M_hyper_mi5_corr_FINAL.rds")

# load: hyper volumes
#HON_H_hyper_mi <- readRDS("HON_H_hyper_mi5.rds")
#HON_M_hyper_mi <- readRDS("HON_M_hyper_mi5.rds")


# set hypervolume comparisons 

set.seed(3)
bm_set_mi_HONHvsHONM <- hypervolume::hypervolume_set(HON_H_hyper_mi, HON_M_hyper_mi, check.memory = FALSE)

# overlap statistics
bm_over_mi <- hypervolume::hypervolume_overlap_statistics(bm_set_mi_HONHvsHONM)

# summarise volumes Historical  vs Modern

hypervolume::get_volume(bm_set_mi_HONHvsHONM)

# Plot 3D hypervolumes 
#Historical vs Modern
bm_set_mi_HONHvsHONM@HVList$HV1 <- NULL
bm_set_mi_HONHvsHONM@HVList$HV2 <- NULL
bm_set_mi_HONHvsHONM@HVList$Union <- NULL

# trait names
names <- c("PC1", "PC2","PC3", "PC4","PC5")

colnames(bm_set_mi_HONHvsHONM@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_HONHvsHONM@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_HONHvsHONM@HVList$Unique_2@RandomPoints) <- names

# plot hypervolumes 
hypervolume::plot.HypervolumeList(bm_set_mi_HONHvsHONM, show.3d = TRUE, plot.3d.axes.id = c(1,2,3), show.random = FALSE, show.data = TRUE, show.centroid = FALSE, show.density = TRUE, cex.random = 5, colors = c("azure3", "red","cornflowerblue"))

#Fancy plot 3D Hypervolume
######################################
#Data frame for Plotly Hyper Volume Plot
inter <- data.frame(bm_set_mi_HONHvsHONM@HVList$Intersection@RandomPoints)
inter <- inter %>% 
  mutate(section = "intersection")

uni1 <- data.frame(bm_set_mi_HONHvsHONM@HVList$Unique_1@RandomPoints) 
uni1 <- uni1 %>% 
  mutate(section = "Unique Historical")

uni2 <- data.frame(bm_set_mi_HONHvsHONM@HVList$Unique_2@RandomPoints)
uni2 <- uni2 %>% 
  mutate(section = "Unique Modern")

HV <- rbind(inter, uni1, uni2)
write.csv(HV,"HyperVolume_HON.csv")

colnames(bm_set_mi_HONHvsHONM@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_HONHvsHONM@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_HONHvsHONM@HVList$Unique_2@RandomPoints) <- names


library(plotly)

HV <- read.csv("HyperVolume_HON.csv", h = T)

pal <- c("azure3","red", "cornflowerblue" )


fig <- plot_ly(HV, x = ~PC1, y = ~(-PC2), z = ~(-PC3), type = 'scatter3d', mode = 'markers', color = ~section, colors = pal, marker = list(size = 3,opacity = 0.6, line = list(color = 'grey66', width = 2)))

fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                   yaxis = list(title = 'PC2'),
                                   zaxis = list(title = 'PC3')))

fig

#FLORENCIA
set.seed(3)
FLOR_H_hyper_mi <- hypervolume::hypervolume_svm(PCA_FLOR_H_A, name = "FLOR_H")

saveRDS(FLOR_H_hyper_mi, "FLOR_H_hyper_mi5_corr_FINAL.rds")

set.seed(3)
FLOR_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_FLOR_M_A, name = "FLOR_M")

saveRDS(FLOR_M_hyper_mi, "FLOR_M_hyper_mi5_corr_FINAL.rds")

# load: hyper volumes
#FLOR_H_hyper_mi <- readRDS("FLOR_H_hyper_mi5.rds")
#FLOR_M_hyper_mi <- readRDS("FLOR_M_hyper_mi5.rds")


# set hypervolume comparisons 

set.seed(3)
bm_set_mi_FLORHvsFLORM <- hypervolume::hypervolume_set(FLOR_H_hyper_mi, FLOR_M_hyper_mi, check.memory = FALSE)

# overlap statistics
bm_over_mi <- hypervolume::hypervolume_overlap_statistics(bm_set_mi_FLORHvsFLORM)

# summarise volumes Historical  vs Modern

hypervolume::get_volume(bm_set_mi_FLORHvsFLORM)

# Plot 3D hypervolumes 
#Historical vs Modern
bm_set_mi_FLORHvsFLORM@HVList$HV1 <- NULL
bm_set_mi_FLORHvsFLORM@HVList$HV2 <- NULL
bm_set_mi_FLORHvsFLORM@HVList$Union <- NULL

# trait names
names <- c("PC1", "PC2","PC3", "PC4","PC5")

colnames(bm_set_mi_FLORHvsFLORM@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_FLORHvsFLORM@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_FLORHvsFLORM@HVList$Unique_2@RandomPoints) <- names

# plot hypervolumes 
hypervolume::plot.HypervolumeList(bm_set_mi_FLORHvsFLORM, show.3d = TRUE, plot.3d.axes.id = c(1,2,3), show.random = FALSE, show.data = TRUE, show.centroid = FALSE, show.density = TRUE, cex.random = 5, colors = c("azure3", "red","cornflowerblue"))

#Fancy plot 3D Hypervolume
######################################
#Data frame for Plotly Hyper Volume Plot
inter <- data.frame(bm_set_mi_FLORHvsFLORM@HVList$Intersection@RandomPoints)
inter <- inter %>% 
  mutate(section = "intersection")

uni1 <- data.frame(bm_set_mi_FLORHvsFLORM@HVList$Unique_1@RandomPoints) 
uni1 <- uni1 %>% 
  mutate(section = "Unique Historical")

uni2 <- data.frame(bm_set_mi_FLORHvsFLORM@HVList$Unique_2@RandomPoints)
uni2 <- uni2 %>% 
  mutate(section = "Unique Modern")

HV <- rbind(inter, uni1, uni2)
write.csv(HV,"HyperVolume_FLOR.csv")

colnames(bm_set_mi_FLORHvsFLORM@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_FLORHvsFLORM@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_FLORHvsFLORM@HVList$Unique_2@RandomPoints) <- names


library(plotly)

HV <- read.csv("HyperVolume_FLOR.csv", h = T)

pal <- c("azure3","red", "cornflowerblue" )


fig <- plot_ly(HV, x = ~PC1, y = ~(-PC2), z = ~(-PC3), type = 'scatter3d', mode = 'markers', color = ~section, colors = pal, marker = list(size = 3,opacity = 0.6, line = list(color = 'grey66', width = 2)))

fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                   yaxis = list(title = 'PC2'),
                                   zaxis = list(title = 'PC3')))

fig

#FUSAGASUGA
set.seed(3)
FUSA_H_hyper_mi <- hypervolume::hypervolume_svm(PCA_FUSA_H_A, name = "FUSA_H")

saveRDS(FUSA_H_hyper_mi, "FUSA_H_hyper_mi5_corr_FINAL.rds")

set.seed(3)
FUSA_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_FUSA_M_A, name = "FUSA_M")

saveRDS(FUSA_M_hyper_mi, "FUSA_M_hyper_mi5_corr_FINAL.rds")

# load: hyper volumes
#FUSA_H_hyper_mi <- readRDS("FUSA_H_hyper_mi5.rds")
#FUSA_M_hyper_mi <- readRDS("FUSA_M_hyper_mi5.rds")


# set hypervolume comparisons 

set.seed(3)
bm_set_mi_FUSAHvsFUSAM <- hypervolume::hypervolume_set(FUSA_H_hyper_mi, FUSA_M_hyper_mi, check.memory = FALSE)

# overlap statistics
bm_over_mi <- hypervolume::hypervolume_overlap_statistics(bm_set_mi_FUSAHvsFUSAM)

# summarise volumes Historical  vs Modern

hypervolume::get_volume(bm_set_mi_FUSAHvsFUSAM)

# Plot 3D hypervolumes 
#Historical vs Modern
bm_set_mi_FUSAHvsFUSAM@HVList$HV1 <- NULL
bm_set_mi_FUSAHvsFUSAM@HVList$HV2 <- NULL
bm_set_mi_FUSAHvsFUSAM@HVList$Union <- NULL

# trait names
names <- c("PC1", "PC2","PC3", "PC4","PC5")

colnames(bm_set_mi_FUSAHvsFUSAM@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_FUSAHvsFUSAM@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_FUSAHvsFUSAM@HVList$Unique_2@RandomPoints) <- names

# plot hypervolumes 
hypervolume::plot.HypervolumeList(bm_set_mi_FUSAHvsFUSAM, show.3d = TRUE, plot.3d.axes.id = c(1,2,3), show.random = FALSE, show.data = TRUE, show.centroid = FALSE, show.density = TRUE, cex.random = 5, colors = c("azure3", "red","cornflowerblue"))

#Fancy plot 3D Hypervolume
######################################
#Data frame for Plotly Hyper Volume Plot
inter <- data.frame(bm_set_mi_FUSAHvsFUSAM@HVList$Intersection@RandomPoints)
inter <- inter %>% 
  mutate(section = "intersection")

uni1 <- data.frame(bm_set_mi_FUSAHvsFUSAM@HVList$Unique_1@RandomPoints) 
uni1 <- uni1 %>% 
  mutate(section = "Unique Historical")

uni2 <- data.frame(bm_set_mi_FUSAHvsFUSAM@HVList$Unique_2@RandomPoints)
uni2 <- uni2 %>% 
  mutate(section = "Unique Modern")

HV <- rbind(inter, uni1, uni2)
write.csv(HV,"HyperVolume_FUSA.csv")

colnames(bm_set_mi_FUSAHvsFUSAM@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_FUSAHvsFUSAM@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_FUSAHvsFUSAM@HVList$Unique_2@RandomPoints) <- names


library(plotly)

HV <- read.csv("HyperVolume_FUSA.csv", h = T)

pal <- c("azure3","red", "cornflowerblue" )


fig <- plot_ly(HV, x = ~PC1, y = ~(-PC2), z = ~(-PC3), type = 'scatter3d', mode = 'markers', color = ~section, colors = pal, marker = list(size = 3,opacity = 0.6, line = list(color = 'grey66', width = 2)))

fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                   yaxis = list(title = 'PC2'),
                                   zaxis = list(title = 'PC3')))

fig
###################################################################


#################################
#2. Functional diversity metrics
####################################################################################
#Functional Diversity Metrics with FD
####################################################################################
library(factoextra)
library(FD)
library(ks)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyverse)
library(data.table)
library(picante)
library(tibble)
library(iNEXT)
library(funrar)

TR <- read.csv("C:\\Users\\camil\\OneDrive - SELVA Investigación para la conservación en el neotropico\\Documentos\\Camila\\Publicaciones\\Chapman_Historical Occupancy\\All_species_traits_CorrJuli_FINAL.csv", sep = ";")

TR <- TR %>% 
  mutate(SD = 1) %>% 
  mutate(log_weight = log(Mass))


#z Transform traits
# function to z-transform data
scale_z <- function(x){
  (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
}

trait_z <- trait %>% 
  mutate(log_weight = log(Mass)) %>% 
  mutate(log_Beak.Length_Culmen = log(Beak.Length_Culmen)) %>% 
  mutate(log_Beak.Width = log(Beak.Width)) %>% 
  mutate(log_Wing.Length = log(Wing.Length)) %>% 
  mutate(log_Tail.Length = log(Tail.Length)) %>% 
  mutate(log_Tarsus.Length = log(Tarsus.Length)) %>% 
  mutate(weight_z = scale_z(log_weight)) %>% 
  mutate(T_cul_z = scale_z(log_Beak.Length_Culmen)) %>% 
  mutate(Bill_w_z = scale_z(log_Beak.Width)) %>% 
  mutate(Wing_l_z = scale_z(log_Wing.Length)) %>% 
  mutate(Tail_l_z = scale_z(log_Tail.Length)) %>% 
  mutate(Tars_z = scale_z(log_Tarsus.Length)) %>% 
  mutate(HWI_z = scale_z(Hand.Wing.Index)) %>% 
  as.data.frame

tr_mi_z <- trait_z %>% 
  select(Species,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z1 <- tr_mi_z[order(tr_mi_z$Species),]

#Species Period Matrix
dat <- read.csv("C:\\Users\\camil\\OneDrive - SELVA Investigación para la conservación en el neotropico\\Documentos\\Camila\\Publicaciones\\Chapman_Historical Occupancy\\SiteSpeciesMatrix_corr.csv", sep = ";",h = T)
head(dat)

#Miraflores

datMIR <- subset(dat,Region == "MIRAFLORES")

datMIR <- datMIR[order(dat$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species, Abundance, Period)

#Remove NAs 
datMIR <- datMIR[complete.cases(datMIR), ]
datMIR <- datMIR[order(datMIR$Species),]

tr_mi_z2 <- datMIR[order(datMIR$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species,Abundance, Period,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z2 <- tr_mi_z2[complete.cases(tr_mi_z2), ]
tr_mi_z2 <- tr_mi_z2[order(tr_mi_z2$Species),] %>% 
  select(Species,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z2 <- tr_mi_z2[!duplicated(tr_mi_z2[,1]), ]

tr_mi_z2 <- data.frame(tr_mi_z2[,-1], row.names=tr_mi_z2[,1])


mat <- datMIR %>% 
  pivot_wider(names_from =  Species, values_from = Abundance) %>% remove_rownames %>%
  column_to_rownames(var = "Period") 

mat[is.na(mat)] <- 0

#Estimating indices of functional diversity

fd_iNDICES <- dbFD(tr_mi_z2, mat,w.abun = FALSE, calc.FRic = TRUE, m = 7, stand.FRic = TRUE , calc.FDiv = TRUE,calc.CWM = FALSE, print.pco = FALSE, messages = TRUE)

fd_iNDICES



#Barbacoas
datBARB <- subset(dat,Region == "BARBACOAS")

datBARB <- datBARB[order(dat$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species, Abundance, Period)

#Remove NAs 
datBARB <- datBARB[complete.cases(datBARB), ]
datBARB <- datBARB[order(datBARB$Species),]

tr_mi_z2 <- datBARB[order(datBARB$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species,Abundance, Period,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z2 <- tr_mi_z2[complete.cases(tr_mi_z2), ]
tr_mi_z2 <- tr_mi_z2[order(tr_mi_z2$Species),] %>% 
  select(Species,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z2 <- tr_mi_z2[!duplicated(tr_mi_z2[,1]), ]

tr_mi_z2 <- data.frame(tr_mi_z2[,-1], row.names=tr_mi_z2[,1])


mat <- datBARB %>% 
  pivot_wider(names_from =  Species, values_from = Abundance) %>% remove_rownames %>%
  column_to_rownames(var = "Period") 

mat[is.na(mat)] <- 0

#Estimating indices of functional diversity

fd_iNDICES <- dbFD(tr_mi_z2, mat,w.abun = FALSE, calc.FRic = TRUE, m = 7, stand.FRic = TRUE , calc.FDiv = TRUE,calc.CWM = FALSE, print.pco = FALSE, messages = TRUE)

fd_iNDICES

#Toche
datTOCHE <- subset(dat,Region == "TOCHE")

datTOCHE <- datTOCHE[order(dat$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species, Abundance, Period)

#Remove NAs 
datTOCHE <- datTOCHE[complete.cases(datTOCHE), ]
datTOCHE <- datTOCHE[order(datTOCHE$Species),]

tr_mi_z2 <- datTOCHE[order(dat$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species,Abundance, Period,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z2 <- tr_mi_z2[complete.cases(tr_mi_z2), ]
tr_mi_z2 <- tr_mi_z2[order(tr_mi_z2$Species),] %>% 
  select(Species,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z2 <- tr_mi_z2[!duplicated(tr_mi_z2[,1]), ]

tr_mi_z2 <- data.frame(tr_mi_z2[,-1], row.names=tr_mi_z2[,1])


mat <- datTOCHE %>% 
  pivot_wider(names_from =  Species, values_from = Abundance) %>% remove_rownames %>%
  column_to_rownames(var = "Period") 

mat[is.na(mat)] <- 0

#Estimating indices of functional diversity

fd_iNDICES <- dbFD(tr_mi_z2, mat,w.abun = FALSE, calc.FRic = TRUE, m = 7, stand.FRic = TRUE , calc.FDiv = TRUE,calc.CWM = FALSE, print.pco = FALSE, messages = TRUE)

fd_iNDICES

#Fusagasuga
datFUSA <- subset(dat,Region == "FUSAGASUGA")

datFUSA <- datFUSA[order(dat$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species, Abundance, Period)

#Remove NAs 
datFUSA <- datFUSA[complete.cases(datFUSA), ]
datFUSA <- datFUSA[order(datFUSA$Species),]

tr_mi_z2 <- datFUSA[order(dat$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species,Abundance, Period,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z2 <- tr_mi_z2[complete.cases(tr_mi_z2), ]
tr_mi_z2 <- tr_mi_z2[order(tr_mi_z2$Species),] %>% 
  select(Species,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z2 <- tr_mi_z2[!duplicated(tr_mi_z2[,1]), ]

tr_mi_z2 <- data.frame(tr_mi_z2[,-1], row.names=tr_mi_z2[,1])


mat <- datFUSA %>% 
  pivot_wider(names_from =  Species, values_from = Abundance) %>% remove_rownames %>%
  column_to_rownames(var = "Period") 

mat[is.na(mat)] <- 0

#Estimating indices of functional diversity

fd_iNDICES <- dbFD(tr_mi_z2, mat,w.abun = FALSE, calc.FRic = TRUE, m = 7, stand.FRic = TRUE , calc.FDiv = TRUE,calc.CWM = FALSE, print.pco = FALSE, messages = TRUE)

fd_iNDICES

#Honda
datHON <- subset(dat,Region == "HONDA")

datHON <- datHON[order(dat$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species, Abundance, Period)

#Remove NAs 
datHON <- datHON[complete.cases(datHON), ]
datHON <- datHON[order(datHON$Species),]

tr_mi_z2 <- datHON[order(dat$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species,Abundance, Period,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z2 <- tr_mi_z2[complete.cases(tr_mi_z2), ]
tr_mi_z2 <- tr_mi_z2[order(tr_mi_z2$Species),] %>% 
  select(Species,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z2 <- tr_mi_z2[!duplicated(tr_mi_z2[,1]), ]

tr_mi_z2 <- data.frame(tr_mi_z2[,-1], row.names=tr_mi_z2[,1])


mat <- datHON %>% 
  pivot_wider(names_from =  Species, values_from = Abundance) %>% remove_rownames %>%
  column_to_rownames(var = "Period") 

mat[is.na(mat)] <- 0

#Estimating indices of functional diversity

fd_iNDICES <- dbFD(tr_mi_z2, mat,w.abun = FALSE, calc.FRic = TRUE, m = 7, stand.FRic = TRUE , calc.FDiv = TRUE,calc.CWM = FALSE, print.pco = FALSE, messages = TRUE)

fd_iNDICES


#San Agustin
datSAGU <- subset(dat,Region == "SAN AGUSTIN")

datSAGU <- datSAGU[order(dat$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species, Abundance, Period)

#Remove NAs 
datSAGU <- datSAGU[complete.cases(datSAGU), ]
datSAGU <- datSAGU[order(datSAGU$Species),]

tr_mi_z2 <- datSAGU[order(dat$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species,Abundance, Period,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z2 <- tr_mi_z2[complete.cases(tr_mi_z2), ]
tr_mi_z2 <- tr_mi_z2[order(tr_mi_z2$Species),] %>% 
  select(Species,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z2 <- tr_mi_z2[!duplicated(tr_mi_z2[,1]), ]

tr_mi_z2 <- data.frame(tr_mi_z2[,-1], row.names=tr_mi_z2[,1])


mat <- datSAGU %>% 
  pivot_wider(names_from =  Species, values_from = Abundance) %>% remove_rownames %>%
  column_to_rownames(var = "Period") 

mat[is.na(mat)] <- 0

#Estimating indices of functional diversity

fd_iNDICES <- dbFD(tr_mi_z2, mat,w.abun = FALSE, calc.FRic = TRUE, m = 7, stand.FRic = TRUE , calc.FDiv = TRUE,calc.CWM = FALSE, print.pco = FALSE, messages = TRUE)

fd_iNDICES

#Florencia
datFLOR <- subset(dat,Region == "FLORENCIA")

datFLOR <- datFLOR[order(dat$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species, Abundance, Period)

#Remove NAs 
datFLOR <- datFLOR[complete.cases(datFLOR), ]
datFLOR <- datFLOR[order(datFLOR$Species),]

tr_mi_z2 <- datFLOR[order(dat$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species,Abundance, Period,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z2 <- tr_mi_z2[complete.cases(tr_mi_z2), ]
tr_mi_z2 <- tr_mi_z2[order(tr_mi_z2$Species),] %>% 
  select(Species,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z2 <- tr_mi_z2[!duplicated(tr_mi_z2[,1]), ]

tr_mi_z2 <- data.frame(tr_mi_z2[,-1], row.names=tr_mi_z2[,1])


mat <- datFLOR %>% 
  pivot_wider(names_from =  Species, values_from = Abundance) %>% remove_rownames %>%
  column_to_rownames(var = "Period") 

mat[is.na(mat)] <- 0

#Estimating indices of functional diversity

fd_iNDICES <- dbFD(tr_mi_z2, mat,w.abun = FALSE, calc.FRic = TRUE, m = 7, stand.FRic = TRUE , calc.FDiv = TRUE,calc.CWM = FALSE, print.pco = FALSE, messages = TRUE)

fd_iNDICES

#3. TPD and species group comparisons (EXTIRPATED VS NEW)
#####################################################################################

library(TPD)

TR <- read.csv("C:\\Users\\camil\\OneDrive - SELVA Investigación para la conservación en el neotropico\\Documentos\\Camila\\Publicaciones\\Chapman_Historical Occupancy\\All_species_traits_CorrJuli_FINAL.csv",sep = ";" ,h = T)

### add scaled trait measures to traits table
trait_z <- TR %>% 
  mutate(log_weight = log(Mass)) %>% 
  mutate(log_Beak.Length_Culmen = log(Beak.Length_Culmen)) %>% 
  mutate(log_Beak.Width = log(Beak.Width)) %>% 
  mutate(log_Wing.Length = log(Wing.Length)) %>% 
  mutate(log_Tail.Length = log(Tail.Length)) %>% 
  mutate(log_Tarsus.Length = log(Tarsus.Length)) %>% 
  mutate(weight_z = scale_z(log_weight)) %>% 
  mutate(T_cul_z = scale_z(log_Beak.Length_Culmen)) %>% 
  mutate(Bill_w_z = scale_z(log_Beak.Width)) %>% 
  mutate(Wing_l_z = scale_z(log_Wing.Length)) %>% 
  mutate(Tail_l_z = scale_z(log_Tail.Length)) %>% 
  mutate(Tars_z = scale_z(log_Tarsus.Length)) %>% 
  mutate(HWI_z = scale_z(Hand.Wing.Index)) %>% 
  as.data.frame

#subset scaled measures
trait_z <- trait_z %>% 
  select(Species,MODERN, HISTORICAL, MIR_H_corr,MIR_M,MIR_EXT,MIR_NEW,POP_MIR,BARB_H_corr,BARB_M,BARB_EXT,BARB_NEW,POP_BARB,SAGU_H_corr, SAGU_M,SAGU_EXT,SAGU_NEW,POP_SAGU, TOCHE_H_corr,TOCHE_M,TOCHE_EXT, TOCHE_NEW,POP_TOCHE, HON_H_corr,HON_M,HON_EXT, HON_NEW,POP_HON, FLOR_H_corr, FLOR_M,FLOR_EXT, FLOR_NEW,POP_FLOR,FUSA_H_corr, FUSA_M,FUSA_EXT, FUSA_NEW, POP_FUSA,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) 

MIR <- trait_z %>% 
  filter(MIR_H_corr >= 0)

BARB <- trait_z %>% 
  filter(BARB_H_corr >= 0)

SAGU <- trait_z %>% 
  filter(SAGU_H_corr >= 0)

TOCHE <- trait_z %>% 
  filter(TOCHE_H_corr >= 0)

HON <- trait_z %>% 
  filter(HON_H_corr >= 0)

FLOR <- trait_z %>% 
  filter(FLOR_H_corr >= 0)

FUSA <- trait_z %>% 
  filter(FUSA_H_corr >= 0)


#### --------------------------------------------------------------
## PCA ## MIRAFLORES
#### --------------------------------------------------------------

#table just with columns of morphological traits and species names as rownames
tr_mi_z_MIR <- MIR %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)%>% 
  column_to_rownames(var = "Species")

#Run PCA
pcaTotal_MIR <- princomp(tr_mi_z_MIR, cor = TRUE, scores = TRUE)
summary(pcaTotal_MIR)

#Scores
scoresTotal_MIR <- as.data.frame(pcaTotal_MIR$scores) %>% 
  tibble::rownames_to_column("Species") %>% 
  left_join(MIR, by = "Species")

#Community TPDc using PCAs as traits 

#Create new database for community TPD to compare extirpated, new additions and share specieswith Total PCA scores
TotalPCA<- scoresTotal_MIR %>% 
  mutate(SD = 1) %>%
  select(Species,Comp.1, Comp.2, Comp.3,POP_MIR,MODERN, HISTORICAL, SD) #POP = groups of species, A = Extirpated, B = Novel additions, C = Shared species

#PC1
TRA <- matrix(c(TotalPCA$Comp.1), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_MIR
ABUN <- TotalPCA %>% 
  select(Species, MODERN,POP_MIR) %>% 
  pivot_wider(names_from = Species, values_from = MODERN) %>% 
  column_to_rownames('POP_MIR')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC2
TRA <- matrix(c(TotalPCA$Comp.2), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_MIR
ABUN <- TotalPCA %>% 
  select(Species, HISTORICAL,POP_MIR) %>% 
  pivot_wider(names_from = Species, values_from = HISTORICAL) %>% 
  column_to_rownames('POP_MIR')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC3

TRA <- matrix(c(TotalPCA$Comp.3), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_MIR
ABUN <- TotalPCA %>% 
  select(Species, HISTORICAL,POP_MIR) %>% 
  pivot_wider(names_from = Species, values_from = HISTORICAL) %>% 
  column_to_rownames('POP_MIR')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN)

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))

#Density comparisons with kolmogorov-smirnov test
library("sm")

A <- subset(TotalPCA, TotalPCA$POP == "A")
B <- subset(TotalPCA, TotalPCA$POP == "B")
C <- subset(TotalPCA, TotalPCA$POP == "C")

#KS test
ks.test(A$Comp.1, B$Comp.1)#P = 0.117

ks.test(A$Comp.2, B$Comp.2)# P = 0.69

ks.test(A$Comp.3, B$Comp.3)# P = 0.04


##Graphs of TPDs for PC1, 2 and 3

TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2) %>% 
  mutate(Comp3M = Comp.3)

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_MIR))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_MIR", values = c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (67%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.8, 0.9))+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "Novel additions", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_MIR))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_MIR", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (20%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "New additions", "Shared species"))

p


#### --------------------------------------------------------------
## PCA ## BARBACOAS
#### --------------------------------------------------------------

#table just with columns of morphological traits and species names as rownames
tr_mi_z_BARB <- BARB %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)%>% 
  column_to_rownames(var = "Species")

#Run PCA
pcaTotal_BARB <- princomp(tr_mi_z_BARB, cor = TRUE, scores = TRUE)
summary(pcaTotal_BARB)

#Scores
scoresTotal_BARB <- as.data.frame(pcaTotal_BARB$scores) %>% 
  tibble::rownames_to_column("Species") %>% 
  left_join(BARB, by = "Species")

#Community TPDc using PCAs as traits 

#Create new database for community TPD to compare extirpated, new additions and share specieswith Total PCA scores
TotalPCA<- scoresTotal_BARB %>% 
  mutate(SD = 1) %>%
  select(Species,Comp.1, Comp.2, Comp.3,POP_BARB,MODERN, HISTORICAL, SD) #POP = groups of species, A = Extirpated, B = Novel additions, C = Shared species

#PC1
TRA <- matrix(c(TotalPCA$Comp.1), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_BARB
ABUN <- TotalPCA %>% 
  select(Species, MODERN,POP_BARB) %>% 
  pivot_wider(names_from = Species, values_from = MODERN) %>% 
  column_to_rownames('POP_BARB')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC2
TRA <- matrix(c(TotalPCA$Comp.2), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_BARB
ABUN <- TotalPCA %>% 
  select(Species, HISTORICAL,POP_BARB) %>% 
  pivot_wider(names_from = Species, values_from = HISTORICAL) %>% 
  column_to_rownames('POP_BARB')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC3

TRA <- matrix(c(TotalPCA$Comp.3), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_BARB
ABUN <- TotalPCA %>% 
  select(Species, HISTORICAL,POP_BARB) %>% 
  pivot_wider(names_from = Species, values_from = HISTORICAL) %>% 
  column_to_rownames('POP_BARB')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN)

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))

#Density comparisons with kolmogorov-smirnov test
library("sm")

A <- subset(TotalPCA, TotalPCA$POP == "A")
B <- subset(TotalPCA, TotalPCA$POP == "B")
C <- subset(TotalPCA, TotalPCA$POP == "C")

#KS test
ks.test(A$Comp.1, B$Comp.1)#P = 0.03

ks.test(A$Comp.2, B$Comp.2)# P = 0.65

ks.test(A$Comp.3, B$Comp.3)# P = 0.34


##Graphs of TPDs for PC1, 2 and 3

TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2) %>% 
  mutate(Comp3M = Comp.3)

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_BARB))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_BARB", values = c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (63%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.8, 0.9))+ geom_text(x = 5, y = 0.2, label = " p = 0.03" , size = 6)+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "Novel additions", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_BARB))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_BARB", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (19%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ geom_text(x = -3, y = 0.4, label = " p = 0.65", size = 6 )+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "New additions", "Shared species"))

p


#### --------------------------------------------------------------
## PCA ## TOCHE
#### --------------------------------------------------------------

#table just with columns of morphological traits and species names as rownames
tr_mi_z_TOCHE <- TOCHE %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)%>% 
  column_to_rownames(var = "Species")

#Run PCA
pcaTotal_TOCHE <- princomp(tr_mi_z_TOCHE, cor = TRUE, scores = TRUE)
summary(pcaTotal_TOCHE)

#Scores
scoresTotal_TOCHE <- as.data.frame(pcaTotal_TOCHE$scores) %>% 
  tibble::rownames_to_column("Species") %>% 
  left_join(TOCHE, by = "Species")

#Community TPDc using PCAs as traits 

#Create new database for community TPD to compare extirpated, new additions and share specieswith Total PCA scores
TotalPCA<- scoresTotal_TOCHE %>% 
  mutate(SD = 1) %>%
  select(Species,Comp.1, Comp.2, Comp.3,POP_TOCHE,MODERN, HISTORICAL, SD) #POP = groups of species, A = Extirpated, B = Novel additions, C = Shared species

#PC1
TRA <- matrix(c(TotalPCA$Comp.1), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_TOCHE
ABUN <- TotalPCA %>% 
  select(Species, MODERN,POP_TOCHE) %>% 
  pivot_wider(names_from = Species, values_from = MODERN) %>% 
  column_to_rownames('POP_TOCHE')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC2
TRA <- matrix(c(TotalPCA$Comp.2), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_TOCHE
ABUN <- TotalPCA %>% 
  select(Species, HISTORICAL,POP_TOCHE) %>% 
  pivot_wider(names_from = Species, values_from = HISTORICAL) %>% 
  column_to_rownames('POP_TOCHE')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC3

TRA <- matrix(c(TotalPCA$Comp.3), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_TOCHE
ABUN <- TotalPCA %>% 
  select(Species, HISTORICAL,POP_TOCHE) %>% 
  pivot_wider(names_from = Species, values_from = HISTORICAL) %>% 
  column_to_rownames('POP_TOCHE')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN)

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))

#Density comparisons with kolmogorov-smirnov test
library("sm")

A <- subset(TotalPCA, TotalPCA$POP == "A")
B <- subset(TotalPCA, TotalPCA$POP == "B")
C <- subset(TotalPCA, TotalPCA$POP == "C")

#KS test
ks.test(A$Comp.1, B$Comp.1)#P = 0.46

ks.test(A$Comp.2, B$Comp.2)# P = 0.95

ks.test(A$Comp.3, B$Comp.3)# P = 0.89


##Graphs of TPDs for PC1, 2 and 3

TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2) %>% 
  mutate(Comp3M = Comp.3)

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_TOCHE))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_TOCHE", values = c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (66%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.75, 0.9))+ 
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "Novel additions", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_TOCHE))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_TOCHE", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (20%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "New additions", "Shared species"))

p


#### --------------------------------------------------------------
## PCA ## HONDA
#### --------------------------------------------------------------

#table just with columns of morphological traits and species names as rownames
tr_mi_z_HON <- HON %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)%>% 
  column_to_rownames(var = "Species")

#Run PCA
pcaTotal_HON <- princomp(tr_mi_z_HON, cor = TRUE, scores = TRUE)
summary(pcaTotal_HON)

#Scores
scoresTotal_HON <- as.data.frame(pcaTotal_HON$scores) %>% 
  tibble::rownames_to_column("Species") %>% 
  left_join(HON, by = "Species")

#Community TPDc using PCAs as traits 

#Create new database for community TPD to compare extirpated, new additions and share specieswith Total PCA scores
TotalPCA<- scoresTotal_HON %>% 
  mutate(SD = 1) %>%
  select(Species,Comp.1, Comp.2, Comp.3,POP_HON,MODERN, HISTORICAL, SD) #POP = groups of species, A = Extirpated, B = Novel additions, C = Shared species

#PC1
TRA <- matrix(c(TotalPCA$Comp.1), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_HON
ABUN <- TotalPCA %>% 
  select(Species, MODERN,POP_HON) %>% 
  pivot_wider(names_from = Species, values_from = MODERN) %>% 
  column_to_rownames('POP_HON')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC2
TRA <- matrix(c(TotalPCA$Comp.2), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_HON
ABUN <- TotalPCA %>% 
  select(Species, HISTORICAL,POP_HON) %>% 
  pivot_wider(names_from = Species, values_from = HISTORICAL) %>% 
  column_to_rownames('POP_HON')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC3

TRA <- matrix(c(TotalPCA$Comp.3), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_HON
ABUN <- TotalPCA %>% 
  select(Species, HISTORICAL,POP_HON) %>% 
  pivot_wider(names_from = Species, values_from = HISTORICAL) %>% 
  column_to_rownames('POP_HON')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN)

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))

#Density comparisons with kolmogorov-smirnov test
library("sm")

A <- subset(TotalPCA, TotalPCA$POP == "A")
B <- subset(TotalPCA, TotalPCA$POP == "B")
C <- subset(TotalPCA, TotalPCA$POP == "C")

#KS test
ks.test(A$Comp.1, B$Comp.1)#P = 0.19

ks.test(A$Comp.2, B$Comp.2)# P = 0.32

ks.test(A$Comp.3, B$Comp.3)# P = 0.55


##Graphs of TPDs for PC1, 2 and 3

TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2 ) %>% 
  mutate(Comp3M = Comp.3 )

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_HON))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_HON", values = c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (64%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.75, 0.9))+ 
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "Novel additions", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_HON))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_HON", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (20%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "New additions", "Shared species"))

p


#### --------------------------------------------------------------
## PCA ## FLORENCIA
#### --------------------------------------------------------------

#table just with columns of morphological traits and species names as rownames
tr_mi_z_FLOR <- FLOR %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)%>% 
  column_to_rownames(var = "Species")

#Run PCA
pcaTotal_FLOR <- princomp(tr_mi_z_FLOR, cor = TRUE, scores = TRUE)
summary(pcaTotal_FLOR)

#Scores
scoresTotal_FLOR <- as.data.frame(pcaTotal_FLOR$scores) %>% 
  tibble::rownames_to_column("Species") %>% 
  left_join(FLOR, by = "Species")

#Community TPDc using PCAs as traits 

#Create new database for community TPD to compare extirpated, new additions and share specieswith Total PCA scores
TotalPCA<- scoresTotal_FLOR %>% 
  mutate(SD = 1) %>%
  select(Species,Comp.1, Comp.2, Comp.3,POP_FLOR,MODERN, HISTORICAL, SD) #POP = groups of species, A = Extirpated, B = Novel additions, C = Shared species

#PC1
TRA <- matrix(c(TotalPCA$Comp.1), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_FLOR
ABUN <- TotalPCA %>% 
  select(Species, MODERN,POP_FLOR) %>% 
  pivot_wider(names_from = Species, values_from = MODERN) %>% 
  column_to_rownames('POP_FLOR')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC2
TRA <- matrix(c(TotalPCA$Comp.2), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_FLOR
ABUN <- TotalPCA %>% 
  select(Species, HISTORICAL,POP_FLOR) %>% 
  pivot_wider(names_from = Species, values_from = HISTORICAL) %>% 
  column_to_rownames('POP_FLOR')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC3

TRA <- matrix(c(TotalPCA$Comp.3), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_FLOR
ABUN <- TotalPCA %>% 
  select(Species, HISTORICAL,POP_FLOR) %>% 
  pivot_wider(names_from = Species, values_from = HISTORICAL) %>% 
  column_to_rownames('POP_FLOR')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN)

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))

#Density comparisons with kolmogorov-smirnov test
library("sm")

A <- subset(TotalPCA, TotalPCA$POP == "A")
B <- subset(TotalPCA, TotalPCA$POP == "B")
C <- subset(TotalPCA, TotalPCA$POP == "C")

#KS test
ks.test(A$Comp.1, B$Comp.1)#P = 0.34

ks.test(A$Comp.2, B$Comp.2)# P = 0.07

ks.test(A$Comp.3, B$Comp.3)# P = 0.01


##Graphs of TPDs for PC1, 2 and 3

TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2 ) %>% 
  mutate(Comp3M = Comp.3 )

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_FLOR))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_FLOR", values = c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (65%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.75, 0.9))+ 
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "Novel additions", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_FLOR))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_FLOR", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (18%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "New additions", "Shared species"))

p


#### --------------------------------------------------------------
## PCA ## SAN AGUSTIN
#### --------------------------------------------------------------

#table just with columns of morphological traits and species names as rownames
tr_mi_z_SAGU <- SAGU %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)%>% 
  column_to_rownames(var = "Species")

#Run PCA
pcaTotal_SAGU <- princomp(tr_mi_z_SAGU, cor = TRUE, scores = TRUE)
summary(pcaTotal_SAGU)

#Scores
scoresTotal_SAGU <- as.data.frame(pcaTotal_SAGU$scores) %>% 
  tibble::rownames_to_column("Species") %>% 
  left_join(SAGU, by = "Species")

#Community TPDc using PCAs as traits 

#Create new database for community TPD to compare extirpated, new additions and share specieswith Total PCA scores
TotalPCA<- scoresTotal_SAGU %>% 
  mutate(SD = 1) %>%
  select(Species,Comp.1, Comp.2, Comp.3,POP_SAGU,MODERN, HISTORICAL, SD) #POP = groups of species, A = Extirpated, B = Novel additions, C = Shared species

#PC1
TRA <- matrix(c(TotalPCA$Comp.1), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_SAGU
ABUN <- TotalPCA %>% 
  select(Species, MODERN,POP_SAGU) %>% 
  pivot_wider(names_from = Species, values_from = MODERN) %>% 
  column_to_rownames('POP_SAGU')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC2
TRA <- matrix(c(TotalPCA$Comp.2), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_SAGU
ABUN <- TotalPCA %>% 
  select(Species, HISTORICAL,POP_SAGU) %>% 
  pivot_wider(names_from = Species, values_from = HISTORICAL) %>% 
  column_to_rownames('POP_SAGU')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC3

TRA <- matrix(c(TotalPCA$Comp.3), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_SAGU
ABUN <- TotalPCA %>% 
  select(Species, HISTORICAL,POP_SAGU) %>% 
  pivot_wider(names_from = Species, values_from = HISTORICAL) %>% 
  column_to_rownames('POP_SAGU')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN)

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))

#Density comparisons with kolmogorov-smirnov test
library("sm")

A <- subset(TotalPCA, TotalPCA$POP == "A")
B <- subset(TotalPCA, TotalPCA$POP == "B")
C <- subset(TotalPCA, TotalPCA$POP == "C")

#KS test
ks.test(A$Comp.1, B$Comp.1)#P = 0.49

ks.test(A$Comp.2, B$Comp.2)# P = 0.61

ks.test(A$Comp.3, B$Comp.3)# P = 0.26


##Graphs of TPDs for PC1, 2 and 3

TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2) %>% 
  mutate(Comp3M = Comp.3 )

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_SAGU))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_SAGU", values = c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (65%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.75, 0.9))+ 
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "Novel additions", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_SAGU))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_SAGU", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (20%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "New additions", "Shared species"))

p


#### --------------------------------------------------------------
## PCA ## FUSAGASUGA
#### --------------------------------------------------------------

#table just with columns of morphological traits and species names as rownames
tr_mi_z_FUSA <- FUSA %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)%>% 
  column_to_rownames(var = "Species")

#Run PCA
pcaTotal_FUSA <- princomp(tr_mi_z_FUSA, cor = TRUE, scores = TRUE)
summary(pcaTotal_FUSA)

#Scores
scoresTotal_FUSA <- as.data.frame(pcaTotal_FUSA$scores) %>% 
  tibble::rownames_to_column("Species") %>% 
  left_join(FUSA, by = "Species")

#Community TPDc using PCAs as traits 

#Create new database for community TPD to compare extirpated, new additions and share specieswith Total PCA scores
TotalPCA<- scoresTotal_FUSA %>% 
  mutate(SD = 1) %>%
  select(Species,Comp.1, Comp.2, Comp.3,POP_FUSA,MODERN, HISTORICAL, SD) #POP = groups of species, A = Extirpated, B = Novel additions, C = Shared species

#PC1
TRA <- matrix(c(TotalPCA$Comp.1), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_FUSA
ABUN <- TotalPCA %>% 
  select(Species, MODERN,POP_FUSA) %>% 
  pivot_wider(names_from = Species, values_from = MODERN) %>% 
  column_to_rownames('POP_FUSA')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC2
TRA <- matrix(c(TotalPCA$Comp.2), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_FUSA
ABUN <- TotalPCA %>% 
  select(Species, HISTORICAL,POP_FUSA) %>% 
  pivot_wider(names_from = Species, values_from = HISTORICAL) %>% 
  column_to_rownames('POP_FUSA')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC3

TRA <- matrix(c(TotalPCA$Comp.3), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_FUSA
ABUN <- TotalPCA %>% 
  select(Species, HISTORICAL,POP_FUSA) %>% 
  pivot_wider(names_from = Species, values_from = HISTORICAL) %>% 
  column_to_rownames('POP_FUSA')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN)

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))

#Density comparisons with kolmogorov-smirnov test
library("sm")

A <- subset(TotalPCA, TotalPCA$POP == "A")
B <- subset(TotalPCA, TotalPCA$POP == "B")
C <- subset(TotalPCA, TotalPCA$POP == "C")

#KS test
ks.test(A$Comp.1, B$Comp.1)#P = 0.54

ks.test(A$Comp.2, B$Comp.2)# P = 0.95

ks.test(A$Comp.3, B$Comp.3)# P = 0.38


##Graphs of TPDs for PC1, 2 and 3

TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2 ) %>% 
  mutate(Comp3M = Comp.3 )

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_FUSA))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_FUSA", values = c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (63%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.75, 0.9))+ 
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "Novel additions", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_FUSA))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_FUSA", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (22%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "New additions", "Shared species"))

p


##########################################
#Mapping vulnerable traits - extirpated species
################################################

#### Load data

# load: trait data
trait <- read.csv("C:\\Users\\camil\\OneDrive - SELVA Investigación para la conservación en el neotropico\\Documentos\\Camila\\Publicaciones\\Chapman_Historical Occupancy\\All_species_traits_CorrJuli_FINAL.csv", sep = ";")

head(trait)
# Exploring trait correlations (verify trait column numbers)

just_trait <- trait[, c(48,49,50,51,52,53,54,55,56,57,58)]
head(just_trait)
res <- cor(just_trait)
round(res, 2)

#z Transform traits
# function to z-transform data
scale_z <- function(x){
  (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
}

### add scaled trait measures to traits table
trait_z <- trait %>% 
  mutate(log_weight = log(Mass)) %>% 
  mutate(log_Beak.Length_Culmen = log(Beak.Length_Culmen)) %>% 
  mutate(log_Beak.Width = log(Beak.Width)) %>% 
  mutate(log_Wing.Length = log(Wing.Length)) %>% 
  mutate(log_Tail.Length = log(Tail.Length)) %>% 
  mutate(log_Tarsus.Length = log(Tarsus.Length)) %>% 
  mutate(weight_z = scale_z(log_weight)) %>% 
  mutate(T_cul_z = scale_z(log_Beak.Length_Culmen)) %>% 
  mutate(Bill_w_z = scale_z(log_Beak.Width)) %>% 
  mutate(Wing_l_z = scale_z(log_Wing.Length)) %>% 
  mutate(Tail_l_z = scale_z(log_Tail.Length)) %>% 
  mutate(Tars_z = scale_z(log_Tarsus.Length)) %>% 
  mutate(HWI_z = scale_z(Hand.Wing.Index)) %>% 
  as.data.frame

#subset scaled measures
trait_z <- trait_z %>% 
  select(Species,MODERN, HISTORICAL,MIR_EXT,BARB_EXT,SAGU_EXT,TOCHE_EXT, HON_EXT, FLOR_EXT, FUSA_EXT, FUSA_NEW,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z, EXTINCT) 

MIR <- trait_z %>% 
  filter(MIR_EXT == 1)

BARB <- trait_z %>% 
  filter(BARB_EXT == 1)

SAGU <- trait_z %>% 
  filter(SAGU_EXT == 1)

TOCHE <- trait_z %>% 
  filter(TOCHE_EXT == 1)

HON <- trait_z %>% 
  filter(HON_EXT == 1)

FLOR <- trait_z %>% 
  filter(FLOR_EXT == 1)

FUSA <- trait_z %>% 
  filter(FUSA_EXT == 1)

TOTAL <- trait_z %>% 
  filter(EXTINCT == 1)




#### --------------------------------------------------------------
## PCA ## TOTAL
#### --------------------------------------------------------------

#table just with columns of morphological traits and species names as rownames
tr_mi_z_TOTAL <- trait_z %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)%>% 
  column_to_rownames(var = "Species")

#Run PCA
pcaTotal_TOTAL <- princomp(tr_mi_z_TOTAL, cor = TRUE, scores = TRUE)
summary(pcaTotal_TOTAL)

#Scores
scoresPCATotal_TOTAL <- as.data.frame(pcaTotal_TOTAL$scores) %>% 
  tibble::rownames_to_column("Species")

scoresPCATotal_TOTAL <- scoresPCATotal_TOTAL %>% 
  # convert long to wide
  tidyr::gather(key, value, -Species) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# loadings
loadingsPCATotal_TOTAL <- as.data.frame(unclass(pcaTotal_TOTAL$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()

loadingsPCATotal_TOTAL <- loadingsPCATotal_TOTAL %>% 
  # convert long to wide
  tidyr::gather(key, value, -trait) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# scalar to adjust arrow size
sc_mi <- 7

loadings_sc_TOTAL <- loadingsPCATotal_TOTAL %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")), funs(.*sc_mi)) %>% 
  # variable names
  mutate(trait = c("Bill width","Hand-Wing Index", "Total culmen", "Tail length", "Tarsus", "Body mass", "Wing length"))

#Subset for species in each Region and period

TOTAL_1 <- trait_z %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  right_join(scoresPCATotal_TOTAL, by = "Species") 

# kernel density estimation for each Region and each period
pc_raw_TOTAL <- TOTAL_1 %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "Species")

# optimal bandwidth estimation
hpi_TOTAL <- Hpi(x = pc_raw_TOTAL)

# kernel density estimation Miraflores_Historical  
est_TOTAL <- kde(x = pc_raw_TOTAL, H = hpi_TOTAL, compute.cont = TRUE)  

den_TOTAL <- list(est_TOTAL$eval.points[[1]], est_TOTAL$eval.points[[2]], est_TOTAL$estimate)
names(den_TOTAL) <- c("x", "y", "z")
dimnames(den_TOTAL$z) <- list(den_TOTAL$x, den_TOTAL$y)
dcc_TOTAL <- melt(den_TOTAL$z)

# source: kernel function
## --------------------------------------------------------------
## Name: 2.1.3-kernel-function.R
## Description: Function to calculate kernel density probablities
## Date: October 2017
## Outputs: Function named 'cl'
## Args:
## df = dataframe of kernel density data
## prob = Probabilty level e.g. 0.95 (95% confidence level)
## --------------------------------------------------------------

cl <- function(df, prob) {
  dx <- diff(df$x[1:2])
  dy <- diff(df$y[1:2])
  sz <- sort(df$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}

# 0.5 probability kernel
cl_50_TOTAL <- cl(df = den_TOTAL, prob = 0.50)
# 0.95 probability kernel
cl_95_TOTAL <- cl(df = den_TOTAL, prob = 0.95)
# 0.99 probability kernel
cl_99_TOTAL <- cl(df = den_TOTAL, prob = 0.99)

# save principal component data

PCA_TOTAL <- trait_z %>% 
    select(Species,MIR_EXT,BARB_EXT,SAGU_EXT,TOCHE_EXT, HON_EXT, FLOR_EXT, FUSA_EXT,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z, EXTINCT) %>% 
  left_join(scoresPCATotal_TOTAL, by = "Species") 

#### PCA plots
library(ggplot2)

# colour palette
col_pal <- colorRampPalette(c("grey30","grey60", "grey80", "white"))(200)


### Folowing lines indentify species that correspond to extinction,
#recolonization or new additions in different periods

extir_MIR <- subset(PCA_TOTAL, MIR_EXT == 1)

extir_BARB <- subset(PCA_TOTAL, BARB_EXT == 1)

extir_SAGU <- subset(PCA_TOTAL, SAGU_EXT == 1)

extir_TOCHE <- subset(PCA_TOTAL, TOCHE_EXT == 1)

extir_HON <- subset(PCA_TOTAL, HON_EXT == 1)

extir_FLOR <- subset(PCA_TOTAL, FLOR_EXT == 1)

extir_FUSA <- subset(PCA_TOTAL, FUSA_EXT == 1)

# plot All Extirpated species

pca_plot_TOTAL <- ggplot(dcc_TOTAL, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_TOTAL, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_MIR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_BARB, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "chocolate1") +
  geom_point(data = extir_SAGU, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "darkgoldenrod1") +
  geom_point(data = extir_TOCHE, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "chartreuse3") +
  geom_point(data = extir_HON, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "chartreuse4") +
  geom_point(data = extir_FLOR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "deepskyblue3") +
  geom_point(data = extir_FUSA, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "deepskyblue4") +
  
  
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_TOTAL, colour = "grey30", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_95_TOTAL, colour = "grey60", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_99_TOTAL, colour = "grey70", linewidth = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_TOTAL, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (63%)", y = "PC2 - Dispersal ability (20%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="All Extirpated Species", color="black", size = 5)

# display plot
windows()
pca_plot_TOTAL

# plot All Extirpated species

pca_plot_TOTAL <- ggplot(dcc_TOTAL, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_TOTAL, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_MIR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_BARB, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_SAGU, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_TOCHE, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_HON, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_FLOR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_FUSA, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  
  
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_TOTAL, colour = "grey30", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_95_TOTAL, colour = "grey60", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_99_TOTAL, colour = "grey70", linewidth = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_TOTAL, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (63%)", y = "PC2 - Dispersal ability (20%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="All Extirpated Species", color="black", size = 5)

# display plot
windows()
pca_plot_TOTAL



# plot Extirpated species - Miraflores

pca_plot_TOTAL <- ggplot(dcc_TOTAL, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_TOTAL, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_MIR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_TOTAL, colour = "grey30", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_95_TOTAL, colour = "grey60", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_99_TOTAL, colour = "grey70", linewidth = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_TOTAL, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (63%)", y = "PC2 - Dispersal ability (20%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Miraflores", color="black", size = 5)

# display plot
windows()
pca_plot_TOTAL


# plot Extirpated species - Barbacoas

pca_plot_TOTAL <- ggplot(dcc_TOTAL, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_TOTAL, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_BARB, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_TOTAL, colour = "grey30", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_95_TOTAL, colour = "grey60", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_99_TOTAL, colour = "grey70", linewidth = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_TOTAL, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (63%)", y = "PC2 - Dispersal ability (20%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Barbacoas", color="black", size = 5)

# display plot
windows()
pca_plot_TOTAL


# plot Extirpated species - San Agustin

pca_plot_TOTAL <- ggplot(dcc_TOTAL, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_TOTAL, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_SAGU, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_TOTAL, colour = "grey30", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_95_TOTAL, colour = "grey60", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_99_TOTAL, colour = "grey70", linewidth = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_TOTAL, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (63%)", y = "PC2 - Dispersal ability (20%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="San Agustín", color="black", size = 5)

# display plot
windows()
pca_plot_TOTAL

# plot Extirpated species - Toche

pca_plot_TOTAL <- ggplot(dcc_TOTAL, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_TOTAL, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_TOCHE, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_TOTAL, colour = "grey30", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_95_TOTAL, colour = "grey60", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_99_TOTAL, colour = "grey70", linewidth = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_TOTAL, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (63%)", y = "PC2 - Dispersal ability (20%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Toche", color="black", size = 5)

# display plot
windows()
pca_plot_TOTAL


# plot Extirpated species - Honda

pca_plot_TOTAL <- ggplot(dcc_TOTAL, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_TOTAL, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_HON, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_TOTAL, colour = "grey30", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_95_TOTAL, colour = "grey60", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_99_TOTAL, colour = "grey70", linewidth = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_TOTAL, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (63%)", y = "PC2 - Dispersal ability (20%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Honda", color="black", size = 5)

# display plot
windows()
pca_plot_TOTAL

# plot Extirpated species - Florencia

pca_plot_TOTAL <- ggplot(dcc_TOTAL, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_TOTAL, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_FLOR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_TOTAL, colour = "grey30", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_95_TOTAL, colour = "grey60", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_99_TOTAL, colour = "grey70", linewidth = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_TOTAL, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (63%)", y = "PC2 - Dispersal ability (20%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Florencia", color="black", size = 5)

# display plot
windows()
pca_plot_TOTAL

# plot Extirpated species - Fusagasuga

pca_plot_TOTAL <- ggplot(dcc_TOTAL, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_TOTAL, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_FUSA, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_TOTAL, colour = "grey30", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_95_TOTAL, colour = "grey60", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_99_TOTAL, colour = "grey70", linewidth = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_TOTAL, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC1 - Size (63%)", y = "PC2 - Dispersal ability (20%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Fusagasugá", color="black", size = 5)

# display plot
windows()
pca_plot_TOTAL
