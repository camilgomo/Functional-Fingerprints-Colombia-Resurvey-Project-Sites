#1. Functional space and hypervolume estimation by region
library(pacman)
pacman::p_load(dplyr, plyr, readr, tbible, FD, ade4, cowplot, mice, reshape2, tidyr, ks, hypervolume, alphallhu, purrr, TTR, plotrix, agricolae, psych)

library(factoextra)
library(ggrepel)
library(tibble)

#### Load data

# load: trait data
trait <- read.csv("C:\\Users\\camil\\OneDrive\\Documents\\Camila\\Publicaciones\\Chapman_Historical Occupancy\\All_species_traits_CorrJuli.csv", sep = ";")

head(trait)
# Exploring trait correlations (verify trait column numbers)

just_trait <- trait[, c(43,44,45,46,47,48,49,50,51,52,53)]
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
  mutate(weight_z = scale_z(log_weight)) %>% 
  mutate(T_cul_z = scale_z(Beak.Length_Culmen)) %>% 
  mutate(Bill_w_z = scale_z(Beak.Width)) %>% 
  mutate(Wing_l_z = scale_z(Wing.Length)) %>% 
  mutate(Tail_l_z = scale_z(Tail.Length)) %>% 
  mutate(Tars_z = scale_z(Tarsus.Length)) %>% 
  mutate(HWI_z = scale_z(Hand.Wing.Index)) %>% 
  as.data.frame

#subset scaled measures
trait_z <- trait_z %>% 
  select(Species,MODERN, HISTORICAL, BARB_H,BARB_M,BARB_EXT,BARB_NEW,SAGU_H, SAGU_M,SAGU_EXT,SAGU_NEW, TOCHE_H,TOCHE_M,TOCHE_EXT, TOCHE_NEW,HON_H,HON_M,HON_EXT, HON_NEW,FLOR_H, FLOR_M,FLOR_EXT, FLOR_NEW, FUSA_H, FUSA_M,FUSA_EXT, FUSA_NEW,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) 

BARB <- trait_z %>% 
  filter(BARB_H >= 0)

SAGU <- trait_z %>% 
  filter(SAGU_H >= 0)

TOCHE <- trait_z %>% 
  filter(TOCHE_H >= 0)

HON <- trait_z %>% 
  filter(HON_H >= 0)

FLOR <- trait_z %>% 
  filter(FLOR_H >= 0)

FUSA <- trait_z %>% 
  filter(FUSA_H >= 0)

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

BARB_H <- trait_z %>% 
  filter(BARB_H >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_BARB, by = "Species") 

BARB_M <- trait_z %>% 
  filter(BARB_M >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_BARB, by = "Species")

SAGU_H <- trait_z %>% 
  filter(SAGU_H >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_SAGU, by = "Species") 

SAGU_M <- trait_z %>% 
  filter(SAGU_M >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_SAGU, by = "Species")

TOCHE_H <- trait_z %>% 
  filter(TOCHE_H >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_TOCHE, by = "Species") 

TOCHE_M <- trait_z %>% 
  filter(TOCHE_M >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_TOCHE, by = "Species")

HON_H <- trait_z %>% 
  filter(HON_H >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_HON, by = "Species") 

HON_M <- trait_z %>% 
  filter(HON_M >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_HON, by = "Species")

FLOR_H <- trait_z %>% 
  filter(FLOR_H >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FLOR, by = "Species") 

FLOR_M <- trait_z %>% 
  filter(FLOR_M >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FLOR, by = "Species")

FUSA_H <- trait_z %>% 
  filter(FUSA_H >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FUSA, by = "Species") 

FUSA_M <- trait_z %>% 
  filter(FUSA_M >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FUSA, by = "Species")

# kernel density estimation for each Region and each period

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

PCA_BARB <- trait_z %>% 
  filter(BARB_H >= 0) %>% 
  select(Species,BARB_H, BARB_M,BARB_EXT,BARB_NEW, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_BARB, by = "Species") 

PCA_SAGU <- trait_z %>% 
  filter(SAGU_H >= 0) %>% 
  select(Species,SAGU_H, SAGU_M,SAGU_EXT,SAGU_NEW, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_SAGU, by = "Species") 

PCA_TOCHE <- trait_z %>% 
  filter(TOCHE_H >= 0) %>% 
  select(Species,TOCHE_H, TOCHE_M,TOCHE_EXT,TOCHE_NEW, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_TOCHE, by = "Species") 

PCA_HON <- trait_z %>% 
  filter(HON_H >= 0) %>% 
  select(Species,HON_H,HON_M,HON_EXT, HON_NEW, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_HON, by = "Species") 

PCA_FLOR <- trait_z %>% 
  filter(FLOR_H >= 0) %>% 
  select(Species,FLOR_H, FLOR_M,FLOR_EXT,FLOR_NEW, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FLOR, by = "Species") 

PCA_FUSA <- trait_z %>% 
  filter(FUSA_H >= 0) %>% 
  select(Species,FUSA_H,FUSA_M,FUSA_EXT,FUSA_NEW, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FUSA, by = "Species") 

write.csv(PCA_BARB, file = "PCA_BARB_corr.csv", row.names = FALSE)
write.csv(PCA_SAGU, file = "PCA_SAGU_corr.csv", row.names = FALSE)
write.csv(PCA_TOCHE, file = "PCA_TOCHE_corr.csv", row.names = FALSE)
write.csv(PCA_HON, file = "PCA_HON_corr.csv", row.names = FALSE)
write.csv(PCA_FLOR, file = "PCA_FLOR_corr.csv", row.names = FALSE)
write.csv(PCA_FUSA, file = "PCA_FUSA_corr.csv", row.names = FALSE)

write.csv(scoresPCATotal_BARB, file = "PCA_Total_BARB_corr.csv")
write.csv(scoresPCATotal_TOCHE, file = "PCA_Total_TOCHE_corr.csv")
write.csv(scoresPCATotal_SAGU, file = "PCA_Total_SAGU_corr.csv")
write.csv(scoresPCATotal_HON, file = "PCA_Total_HON_corr.csv")
write.csv(scoresPCATotal_FLOR, file = "PCA_Total_FLOR_corr.csv")
write.csv(scoresPCATotal_FUSA, file = "PCA_Total_FUSA_corr.csv")

write.csv(loadingsPCATotal_BARB, file = "Loadings_Total_BARB_corr.csv", row.names = FALSE)
write.csv(loadingsPCATotal_SAGU, file = "Loadings_Total_SAGU.csv_corr", row.names = FALSE)
write.csv(loadingsPCATotal_TOCHE, file = "Loadings_Total_TOCHE_corr.csv", row.names = FALSE)
write.csv(loadingsPCATotal_HON, file = "Loadings_Total_HON_corr.csv", row.names = FALSE)
write.csv(loadingsPCATotal_FLOR, file = "Loadings_Total_FLOR_corr.csv", row.names = FALSE)
write.csv(loadingsPCATotal_FUSA, file = "Loadings_Total_FUSA_corr.csv", row.names = FALSE)

#### PCA plots
library(ggplot2)

# colour palette
col_pal <- colorRampPalette(c("grey30","grey60", "grey80", "white"))(200)


### Folowing lines indentify species that correspond to extinction,
#recolonization or new additions in different periods

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
  labs(x = "PC1 - Size (61%)", y = "PC2 - Dispersal Hability (16%)") +
  xlim(-5,15) +
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
  labs(x = "PC1 - Size (61%)", y = "PC2 - Dispersal Hability (16%)") +
  xlim(-5,15) +
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
  labs(x = "PC1 - Size (62%)", y = "PC2 - Dispersal Hability (17%)") +
  xlim(-5,15) +
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
  labs(x = "PC1 - Size (62%)", y = "PC2 - Dispersal Hability (17%)") +
  xlim(-5,15) +
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
  labs(x = "PC1 - Size (62%)", y = "PC2 - Dispersal Hability (17%)") +
  xlim(-5,15) +
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
  labs(x = "PC1 - Size (62%)", y = "PC2 - Dispersal Hability (17%)") +
  xlim(-5,15) +
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
  labs(x = "PC1 - Size (62%)", y = "PC2 - Dispersal Hability (16%)") +
  xlim(-5,15) +
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
  labs(x = "PC1 - Size (62%)", y = "PC2 - Dispersal Hability (16%)") +
  xlim(-5,15) +
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
  labs(x = "PC1 - Size (61%)", y = "PC2 - Dispersal Hability (15%)") +
  xlim(-5,15) +
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
  labs(x = "PC1 - Size (61%)", y = "PC2 - Dispersal Hability (15%)") +
  xlim(-5,15) +
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
  labs(x = "PC1 - Size (60%)", y = "PC2 - Dispersal Hability (18%)") +
  xlim(-5,15) +
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
  labs(x = "PC1 - Size (60%)", y = "PC2 - Dispersal Hability (18%)") +
  xlim(-5,15) +
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
#BARBACOAS
PCA_BARB_H_A <- subset(PCA_BARB, BARB_H >= 1)
  
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
PCA_TOCHE_H_A <- subset(PCA_TOCHE, TOCHE_H >= 1)

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
PCA_SAGU_H_A <- subset(PCA_SAGU, SAGU_H >= 1)

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
PCA_HON_H_A <- subset(PCA_HON, HON_H >= 1)

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
PCA_FLOR_H_A <- subset(PCA_FLOR, FLOR_H >= 1)

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
PCA_FUSA_H_A <- subset(PCA_FUSA, FUSA_H >= 1)

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
#Barbacoas 
set.seed(3)
BARB_H_hyper_mi <- hypervolume::hypervolume_svm(PCA_BARB_H_A, name = "BARB_H")

saveRDS(BARB_H_hyper_mi, "BARB_H_hyper_mi5_corr.rds")

set.seed(3)
BARB_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_BARB_M_A, name = "BARB_M")

saveRDS(BARB_M_hyper_mi, "BARB_M_hyper_mi5_corr.rds")

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

saveRDS(TOCHE_H_hyper_mi, "TOCHE_H_hyper_mi5_corr.rds")

set.seed(3)
TOCHE_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_TOCHE_M_A, name = "TOCHE_M")

saveRDS(TOCHE_M_hyper_mi, "TOCHE_M_hyper_mi5_corr.rds")

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

saveRDS(SAGU_H_hyper_mi, "SAGU_H_hyper_mi5_corr.rds")

set.seed(3)
SAGU_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_SAGU_M_A, name = "SAGU_M")

saveRDS(SAGU_M_hyper_mi, "SAGU_M_hyper_mi5_corr.rds")

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

saveRDS(HON_H_hyper_mi, "HON_H_hyper_mi5_corr.rds")

set.seed(3)
HON_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_HON_M_A, name = "HON_M")

saveRDS(HON_M_hyper_mi, "HON_M_hyper_mi5_corr.rds")

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

saveRDS(FLOR_H_hyper_mi, "FLOR_H_hyper_mi5_corr.rds")

set.seed(3)
FLOR_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_FLOR_M_A, name = "FLOR_M")

saveRDS(FLOR_M_hyper_mi, "FLOR_M_hyper_mi5_corr.rds")

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

saveRDS(FUSA_H_hyper_mi, "FUSA_H_hyper_mi5_corr.rds")

set.seed(3)
FUSA_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_FUSA_M_A, name = "FUSA_M")

saveRDS(FUSA_M_hyper_mi, "FUSA_M_hyper_mi5_corr.rds")

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

TR <- read.csv("C:\\Users\\camil\\OneDrive\\Documents\\Camila\\Publicaciones\\Chapman_Historical Occupancy\\All_species_traits_CorrJuli.csv", sep = ";")

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
  mutate(weight_z = scale_z(log_weight)) %>% 
  mutate(T_cul_z = scale_z(Beak.Length_Culmen)) %>% 
  mutate(Bill_w_z = scale_z(Beak.Width)) %>% 
  mutate(Wing_l_z = scale_z(Wing.Length)) %>% 
  mutate(Tail_l_z = scale_z(Tail.Length)) %>% 
  mutate(Tars_z = scale_z(Tarsus.Length)) %>% 
  mutate(HWI_z = scale_z(Hand.Wing.Index)) %>% 
  as.data.frame

tr_mi_z <- trait_z %>% 
  select(Species,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z1 <- tr_mi_z[order(tr_mi_z$Species),]

#Species Period Matrix
dat <- read.csv("C:\\Users\\camil\\OneDrive\\Documents\\Camila\\Publicaciones\\Chapman_Historical Occupancy\\SiteSpeciesMatrix_corr.csv", sep = ";",h = T)
head(dat)

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

TR <- read.csv("C:\\Users\\camil\\OneDrive\\Documents\\Camila\\Publicaciones\\Chapman_Historical Occupancy\\All_species_traits_CorrJuli.csv",sep = ";" ,h = T)

### add scaled trait measures to traits table
trait_z <- TR %>% 
  mutate(log_weight = log(Mass)) %>% 
  mutate(weight_z = scale_z(log_weight)) %>% 
  mutate(T_cul_z = scale_z(Beak.Length_Culmen)) %>% 
  mutate(Bill_w_z = scale_z(Beak.Width)) %>% 
  mutate(Wing_l_z = scale_z(Wing.Length)) %>% 
  mutate(Tail_l_z = scale_z(Tail.Length)) %>% 
  mutate(Tars_z = scale_z(Tarsus.Length)) %>% 
  mutate(HWI_z = scale_z(Hand.Wing.Index)) %>% 
  as.data.frame

#subset scaled measures
trait_z <- trait_z %>% 
  select(Species,MODERN, HISTORICAL, BARB_H,BARB_M,BARB_EXT,BARB_NEW,POP_BARB,SAGU_H, SAGU_M,SAGU_EXT,SAGU_NEW,POP_SAGU, TOCHE_H,TOCHE_M,TOCHE_EXT, TOCHE_NEW,POP_TOCHE, HON_H,HON_M,HON_EXT, HON_NEW,POP_HON, FLOR_H, FLOR_M,FLOR_EXT, FLOR_NEW,POP_FLOR,FUSA_H, FUSA_M,FUSA_EXT, FUSA_NEW, POP_FUSA,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) 

BARB <- trait_z %>% 
  filter(BARB_H >= 0)

SAGU <- trait_z %>% 
  filter(SAGU_H >= 0)

TOCHE <- trait_z %>% 
  filter(TOCHE_H >= 0)

HON <- trait_z %>% 
  filter(HON_H >= 0)

FLOR <- trait_z %>% 
  filter(FLOR_H >= 0)

FUSA <- trait_z %>% 
  filter(FUSA_H >= 0)

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
ks.test(A$Comp.1, B$Comp.1)#P = 0.012

ks.test(A$Comp.2, B$Comp.2)# P = 0.31

ks.test(A$Comp.3, B$Comp.3)# P = 0.88


##Graphs of TPDs for PC1, 2 and 3

#Multiply PC2 and 3 by -1 to ease interpretation of increasing values.
TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2 *-1) %>% 
  mutate(Comp3M = Comp.3 *-1)

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_BARB))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_BARB", values = c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (61%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.75, 0.9))+ geom_text(x = 6, y = 0.2, label = " p = 0.012" , size = 6)+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "Novel additions", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_BARB))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_BARB", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (16%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ geom_text(x = 2, y = 0.6, label = " p = 0.31", size = 6 )+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "New additions", "Shared species"))

p

#PC3
windows()
p = ggplot(TotalPCA1, aes(x = Comp3M,fill = POP_BARB))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_BARB", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC3 (12%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none") + geom_text(x = 5, y = 1.5, label = " p = 0.88", size = 6 )+
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
ks.test(A$Comp.1, B$Comp.1)#P = 0.24

ks.test(A$Comp.2, B$Comp.2)# P = 0.71

ks.test(A$Comp.3, B$Comp.3)# P = 0.15


##Graphs of TPDs for PC1, 2 and 3

#Multiply PC2 and 3 by -1 to ease interpretation of increasing values.
TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2 *-1) %>% 
  mutate(Comp3M = Comp.3 *-1)

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_TOCHE))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_TOCHE", values = c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (62%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.75, 0.9))+ geom_text(x = 6, y = 0.2, label = " p = 0.24" , size = 6)+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "Novel additions", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_TOCHE))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_TOCHE", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (16%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ geom_text(x = 3, y = 0.6, label = " p = 0.71", size = 6 )+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "New additions", "Shared species"))

p

#PC3
windows()
p = ggplot(TotalPCA1, aes(x = Comp3M,fill = POP_TOCHE))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_TOCHE", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC3 (9%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none") + geom_text(x = 3.5, y = 1.5, label = " p = 0.15", size = 6 )+
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
ks.test(A$Comp.1, B$Comp.1)#P = 0.15

ks.test(A$Comp.2, B$Comp.2)# P = 0.16

ks.test(A$Comp.3, B$Comp.3)# P = 0.64


##Graphs of TPDs for PC1, 2 and 3

#Multiply PC2 and 3 by -1 to ease interpretation of increasing values.
TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2 *-1) %>% 
  mutate(Comp3M = Comp.3 *-1)

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_HON))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_HON", values = c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (62%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.75, 0.9))+ geom_text(x = 10, y = 0.2, label = " p = 0.15" , size = 6)+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "Novel additions", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_HON))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_HON", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (16%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ geom_text(x = 3, y = 0.6, label = " p = 0.16", size = 6 )+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "New additions", "Shared species"))

p

#PC3
windows()
p = ggplot(TotalPCA1, aes(x = Comp3M,fill = POP_HON))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_HON", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC3 (9%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none") + geom_text(x = 4, y = 0.9, label = " p = 0.64", size = 6 )+
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
ks.test(A$Comp.1, B$Comp.1)#P = 0.29

ks.test(A$Comp.2, B$Comp.2)# P = 0.02

ks.test(A$Comp.3, B$Comp.3)# P = 0.03


##Graphs of TPDs for PC1, 2 and 3

#Multiply PC2 and 3 by -1 to ease interpretation of increasing values.
TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2 *-1) %>% 
  mutate(Comp3M = Comp.3 *-1)

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_FLOR))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_FLOR", values = c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (61%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.75, 0.9))+ geom_text(x = 6, y = 0.2, label = " p = 0.29" , size = 6)+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "Novel additions", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_FLOR))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_FLOR", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (15%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ geom_text(x = 3, y = 0.7, label = " p = 0.03", size = 6 )+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "New additions", "Shared species"))

p

#PC3
windows()
p = ggplot(TotalPCA1, aes(x = Comp3M,fill = POP_FLOR))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_FLOR", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC3 (10%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none") + geom_text(x = 5, y = 1.25, label = " p = 0.03", size = 6 )+
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
ks.test(A$Comp.1, B$Comp.1)#P = 0.69

ks.test(A$Comp.2, B$Comp.2)# P = 0.40

ks.test(A$Comp.3, B$Comp.3)# P = 0.21


##Graphs of TPDs for PC1, 2 and 3

#Multiply PC2 and 3 by -1 to ease interpretation of increasing values.
TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2 *-1) %>% 
  mutate(Comp3M = Comp.3 *-1)

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_SAGU))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_SAGU", values = c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (62%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.75, 0.9))+ geom_text(x = 6, y = 0.2, label = " p = 0.69" , size = 6)+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "Novel additions", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_SAGU))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_SAGU", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (17%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ geom_text(x = 2.5, y = 0.6, label = " p = 0.40", size = 6 )+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "New additions", "Shared species"))

p

#PC3
windows()
p = ggplot(TotalPCA1, aes(x = Comp3M,fill = POP_SAGU))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_SAGU", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC3 (10%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none") + geom_text(x = 6, y = 1.5, label = " p = 0.21", size = 6 )+
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
ks.test(A$Comp.1, B$Comp.1)#P = 0.60

ks.test(A$Comp.2, B$Comp.2)# P = 0.75

ks.test(A$Comp.3, B$Comp.3)# P = 0.25


##Graphs of TPDs for PC1, 2 and 3

#Multiply PC2 and 3 by -1 to ease interpretation of increasing values.
TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2 *-1) %>% 
  mutate(Comp3M = Comp.3 *-1)

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_FUSA))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_FUSA", values = c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (60%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.75, 0.9))+ geom_text(x = 9, y = 0.25, label = " p = 0.60" , size = 6)+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "Novel additions", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_FUSA))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_FUSA", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (18%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ geom_text(x = 3.5, y = 0.7, label = " p = 0.75", size = 6 )+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "New additions", "Shared species"))

p

#PC3
windows()
p = ggplot(TotalPCA1, aes(x = Comp3M,fill = POP_FUSA))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_FUSA", values=c("red", "cornflowerblue", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC3 (9%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none") + geom_text(x = 6, y = 1.5, label = " p = 0.25", size = 6 )+
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), 
                    name=NULL,breaks=c("A", "B", "C"),labels=c("Extirpated species", "New additions", "Shared species"))

p


#########################

#4. Functional Uniqueness and distinctiveness
#####################################################################################

library(funrar)

# load data

#Species Period Matrix
dat <- read.csv("C:\\Users\\camil\\OneDrive\\Documents\\Camila\\Publicaciones\\Chapman_Historical Occupancy\\SiteSpeciesMatrix_corr.csv", sep = ";", h = T)

head(dat)

#Barbacoas
datBARB <- subset(dat,Region == "BARBACOAS")

datBARB <- datBARB[order(dat$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species, Abundance, Period)

#Remove NAs 
datBARB <- datBARB[complete.cases(datBARB), ]
datBARB <- datBARB[order(datBARB$Species),]

tr_mi_z2 <- datBARB[order(dat$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species,Abundance, Period,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z2 <- tr_mi_z2[complete.cases(tr_mi_z2), ]
tr_mi_z2 <- tr_mi_z2[order(tr_mi_z2$Species),]

tr_mi_z2 <- tr_mi_z2[!duplicated(tr_mi_z2[,1]), ]

tr_mi_z2 <- data.frame(tr_mi_z2[,-1], row.names=tr_mi_z2[,1])


matBARB <- datBARB %>% 
  pivot_wider(names_from =  Species, values_from = Abundance) %>% remove_rownames %>%
  column_to_rownames(var = "Period") 

matBARB[is.na(matBARB)] <- 0

#BARBACOAS ###############

matBARB <- as.matrix(matBARB)

#Compute distance matrix
dist_mat <- compute_dist_matrix(tr_mi_z_BARB, metric = "euclidean")

# Compute distinctiveness for each species on each site
di_df = distinctiveness_stack(com_df = datBARB,  # The site x species table
                              sp_col = "Species",  # Name of the species column
                              com = "Period",  # Name of the community column
                              #abund = "Abundance",  # Relative abundances column (facultative)
                              dist_matrix = dist_mat)  # Functional Distance matrix

head(di_df)
write.csv(di_df, file = "Distinctiveness_BARB_corr.csv", row.names = FALSE)

Di_Sum <- di_df %>% 
  group_by(Species) %>%
  summarise_at(vars(Di), funs(mean(., na.rm=TRUE)))

Di_Sum <- Di_Sum[order(Di_Sum$Di),] %>% 
  mutate(ID = c(1:208)) %>% 
  left_join(trait, by = "Species")

Di_Sum <- Di_Sum %>% 
  filter(!is.na(POP_BARB))
########
#Randomizartion to compare means between extirpated and new

#Observed means
A <- subset(Di_Sum, Di_Sum$POP_BARB == "A")
B <- subset(Di_Sum, Di_Sum$POP_BARB == "B")
C <- subset(Di_Sum, Di_Sum$POP_BARB == "C")
meanA <- mean(A$Di)
meanA #2.60

meanB <- mean(B$Di)
meanB #3.29

meanC <- mean(C$Di)
meanC #2.69

obsDifAB <- meanA - meanB #-0.6935
obsDifBC <- meanB - meanC #0.6048


iter    <- 999;          # Total iterations (+1 for observed data = 10k)
diff_AB    <- NULL;# To add difference between groups
diff_BC    <- NULL;# To add difference between groups

N_obs <- 208; # Total number of observations


for(i in 1:iter){   
  reshufled <- Di_Sum
  reshufled$POP_BARB   <- sample(reshufled$POP_BARB, N_obs, replace = FALSE);
  
  mean_A <- mean(reshufled %>% filter(POP_BARB == "A") %>% pull(Di))
  mean_B <- mean(reshufled %>% filter(POP_BARB == "B") %>% pull(Di))
  mean_C <- mean(reshufled %>% filter(POP_BARB == "C") %>% pull(Di))
  
  diff_simAB <- mean_A - mean_B
  diff_simBC <- mean_B - mean_C
  
  diff_AB[i] <- diff_simAB
  diff_BC[i] <- diff_simBC
}

#P value Extirpated  -  New
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_AB), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifAB, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_AB <= obsDifAB) + 1;
total_generated   <- length(diff_AB) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.003

#P value New - Shared
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_BC), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifBC, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_BC >= obsDifBC) + 1;
total_generated   <- length(diff_BC) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.003

### GRaph POP
windows()
p <- ggplot(Di_Sum, aes(x=ID.x, y=Di, fill = POP_BARB, na.rm = T)) + 
  geom_bar(stat = "identity")+ 
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), breaks=c("A","B","C"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional distinctiveness - Di")

p
#Boxplot  Di POP

windows()
p <- ggplot(Di_Sum, aes(x=POP_BARB, y=Di, fill = POP_BARB, na.rm = T)) + 
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(shape=16,position=position_jitter(0.2), colour = "black") +
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), labels = c("Extirpated","New additions","Shared species"))+
  scale_x_discrete(labels = c("Extirpated","New additions","Shared species"))+
  ylim(0,12)+
  #Text for significant differences and p values between group. Replace values with comparison results.
  geom_text(x = 1.5, y = 11.5, label = "p = 0.002")+
  geom_text(x = 2.5, y = 11.5, label = "p = 0.003")+
  
  geom_segment(aes(x = 1, y = 11, xend = 1.9, yend = 11))+
  geom_segment(aes(x = 2, y = 11, xend = 2.9, yend = 11))+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional Distinctiveness")

p

#######GRaph F Distinctiveness by period


Y1 <- subset(di_df, di_df$Period == "HISTORICAL")

Y2 <- subset(di_df, di_df$Period == "MODERN")

meanY1 <- mean(Y1$Di)
meanY2 <- mean(Y2$Di)

d <- di_df %>% 
  group_by(Period) %>%
  summarise_at(vars(Di), funs(mean(., na.rm=TRUE)))

windows()
p <- ggplot(di_df, aes(x=as.factor(Period), y=Di)) + 
  geom_violin(trim = FALSE, colour = "gray40")+ 
  geom_jitter(shape=16, position=position_jitter(0.2), colour = "gray40") +
  geom_jitter(data = Y1, aes(x = as.factor(Period), y = Di), size = 3, alpha = 0.5, colour = "red",position=position_jitter(0.2)) + 
  geom_jitter(data = Y2, aes(x = as.factor(Period), y = Di), size = 3, alpha = 0.5, colour = "cornflowerblue",position=position_jitter(0.2)) + 
  geom_segment(aes(x = 0.5, y = 2.58, xend = 1.5, yend = 2.58))+
  geom_segment(aes(x = 1.5, y = 3.02, xend = 2.5, yend = 3.02))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional distinctiveness")+ scale_x_discrete(labels = c("Historical", "Modern"))

p


#Compute Uniqueness for each species
ui_df = uniqueness_stack(datBARB, "Species", dist_mat)

head(ui_df)
write.csv(ui_df, file = "Uniqueness_BARB_corr.csv", row.names = FALSE)

### GRaph POP

Ui_Sum <- ui_df %>% 
  group_by(Species) %>%
  summarise_at(vars(Ui), funs(mean(., na.rm=TRUE)))

Ui_Sum <- Ui_Sum[order(Ui_Sum$Ui),] %>% 
  mutate(ID = c(1:208)) %>% 
  left_join(trait, by = "Species")

Ui_Sum <- Ui_Sum %>% 
  filter(!is.na(POP_BARB))
########
#Randomizartion to compare means between extirpated and new

#Observed means
A <- subset(Ui_Sum, Ui_Sum$POP_BARB == "A")
B <- subset(Ui_Sum, Ui_Sum$POP_BARB == "B")
C <- subset(Ui_Sum, Ui_Sum$POP_BARB == "C")
meanA <- mean(A$Ui)
meanA #0.6321

meanB <- mean(B$Ui)
meanB #0.6603

meanC <- mean(C$Ui)
meanC #0.4712

obsDifAB <- meanA - meanB #-0.02821
obsDifBC <- meanB - meanC #0.1891

iter    <- 999;          # Total iterations (+1 for observed data = 10k)
diff_AB    <- NULL;# To add difference between groups
diff_BC    <- NULL;# To add difference between groups

N_obs <- 208; # Total number of observations


for(i in 1:iter){   
  reshufled <- Ui_Sum
  reshufled$POP_BARB   <- sample(reshufled$POP_BARB, N_obs, replace = FALSE);
  
  mean_A <- mean(reshufled %>% filter(POP_BARB == "A") %>% pull(Ui))
  mean_B <- mean(reshufled %>% filter(POP_BARB == "B") %>% pull(Ui))
  mean_C <- mean(reshufled %>% filter(POP_BARB == "C") %>% pull(Ui))
  
  diff_simAB <- mean_A - mean_B
  diff_simBC <- mean_B - mean_C
  
  diff_AB[i] <- diff_simAB
  diff_BC[i] <- diff_simBC
}

#P value Extirpated  -  New
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_AB), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifAB, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_AB >= obsDifAB) + 1;
total_generated   <- length(diff_AB) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.61

#P value New - Shared
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_BC), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifBC, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_BC >= obsDifBC) + 1;
total_generated   <- length(diff_BC) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.009


#Boxplot  Ui POP

windows()
p <- ggplot(Ui_Sum, aes(x=POP_BARB, y=Ui, fill = POP_BARB)) + 
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(shape=16,position=position_jitter(0.2), colour = "black") +
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), labels = c("Extirpated","New additions","Shared species"))+
  ylim(0,3)+
  #Text for significant differences and p values between group. Replace values with comparison results.
  geom_text(x = 1.5, y = 3, label = "p = 0.61")+
  geom_text(x = 2.5, y = 3, label = "p = 0.009")+
  
  geom_segment(aes(x = 1, y = 2.7, xend = 1.9, yend = 2.7))+
  geom_segment(aes(x = 2, y = 2.7, xend = 2.9, yend = 2.7))+
  scale_x_discrete(labels = c("Extirpated","New additions","Shared species"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional Uniqueness")

p

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
tr_mi_z2 <- tr_mi_z2[order(tr_mi_z2$Species),]

tr_mi_z2 <- tr_mi_z2[!duplicated(tr_mi_z2[,1]), ]

tr_mi_z2 <- data.frame(tr_mi_z2[,-1], row.names=tr_mi_z2[,1])


matTOCHE <- datTOCHE %>% 
  pivot_wider(names_from =  Species, values_from = Abundance) %>% remove_rownames %>%
  column_to_rownames(var = "Period") 

matTOCHE[is.na(matTOCHE)] <- 0
##################
#TOCHE
matTOCHE <- as.matrix(matTOCHE)


#Compute distance matrix
dist_mat <- compute_dist_matrix(tr_mi_z_TOCHE, metric = "euclidean")

# Compute distinctiveness for each species on each site
di_df = distinctiveness_stack(com_df = datTOCHE,  # The site x species table
                              sp_col = "Species",  # Name of the species column
                              com = "Period",  # Name of the community column
                              #abund = "Abundance",  # Relative abundances column (facultative)
                              dist_matrix = dist_mat)  # Functional Distance matrix

head(di_df)
write.csv(di_df, file = "Distinctiveness_TOCHE_corr.csv", row.names = FALSE)

Di_Sum <- di_df %>% 
  group_by(Species) %>%
  summarise_at(vars(Di), funs(mean(., na.rm=TRUE)))

Di_Sum <- Di_Sum[order(Di_Sum$Di),] %>% 
  mutate(ID = c(1:239)) %>% 
  left_join(trait, by = "Species")

Di_Sum <- Di_Sum %>% 
  filter(!is.na(POP_TOCHE))
########
#Randomizartion to compare means between extirpated and new

#Observed means
A <- subset(Di_Sum, Di_Sum$POP_TOCHE == "A")
B <- subset(Di_Sum, Di_Sum$POP_TOCHE == "B")
C <- subset(Di_Sum, Di_Sum$POP_TOCHE == "C")
meanA <- mean(A$Di)
meanA #2.906

meanB <- mean(B$Di)
meanB #3.197

meanC <- mean(C$Di)
meanC #2.861

obsDifAB <- meanA - meanB #-0.2902
obsDifBC <- meanB - meanC #0.335


iter    <- 999;          # Total iterations (+1 for observed data = 10k)
diff_AB    <- NULL;# To add difference between groups
diff_BC    <- NULL;# To add difference between groups

N_obs <- 239; # Total number of observations


for(i in 1:iter){   
  reshufled <- Di_Sum
  reshufled$POP_TOCHE   <- sample(reshufled$POP_TOCHE, N_obs, replace = FALSE);
  
  mean_A <- mean(reshufled %>% filter(POP_TOCHE == "A") %>% pull(Di))
  mean_B <- mean(reshufled %>% filter(POP_TOCHE == "B") %>% pull(Di))
  mean_C <- mean(reshufled %>% filter(POP_TOCHE == "C") %>% pull(Di))
  
  diff_simAB <- mean_A - mean_B
  diff_simBC <- mean_B - mean_C
  
  diff_AB[i] <- diff_simAB
  diff_BC[i] <- diff_simBC
}

#P value Extirpated  -  New
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_AB), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifAB, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_AB <= obsDifAB) + 1;
total_generated   <- length(diff_AB) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.15

#P value New - Shared
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_BC), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifBC, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_BC >= obsDifBC) + 1;
total_generated   <- length(diff_BC) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.027

### GRaph POP
windows()
p <- ggplot(Di_Sum, aes(x=ID.x, y=Di, fill = POP_TOCHE)) + 
  geom_bar(stat = "identity")+ 
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), breaks=c("A","B","C"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional distinctiveness - Di")

p
#Boxplot  Di POP

windows()
p <- ggplot(Di_Sum, aes(x=POP_TOCHE, y=Di, fill = POP_TOCHE)) + 
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(shape=16,position=position_jitter(0.2), colour = "black") +
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), labels = c("Extirpated","New additions","Shared species"))+
  scale_x_discrete(labels = c("Extirpated","New additions","Shared species"))+
  ylim(0,12)+
  #Text for significant differences and p values between group. Replace values with comparison results.
  geom_text(x = 1.5, y = 11, label = "p = 0.15")+
  geom_text(x = 2.5, y = 11, label = "p = 0.02")+
  
  geom_segment(aes(x = 1, y = 10.5, xend = 1.9, yend = 10.5))+
  geom_segment(aes(x = 2, y = 10.5, xend = 2.9, yend = 10.5))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional Distinctiveness")

p

#######GRaph F Distinctiveness by period

Y1 <- subset(di_df, di_df$Period == "HISTORICAL")

Y2 <- subset(di_df, di_df$Period == "MODERN")

meanY1 <- mean(Y1$Di)
meanY2 <- mean(Y2$Di)

d <- di_df %>% 
  group_by(Period) %>%
  summarise_at(vars(Di), funs(mean(., na.rm=TRUE)))

windows()
p <- ggplot(di_df, aes(x=as.factor(Period), y=Di)) + 
  geom_violin(trim = FALSE, colour = "gray40")+ 
  geom_jitter(shape=16, position=position_jitter(0.2), colour = "gray40") +
  geom_jitter(data = Y1, aes(x = as.factor(Period), y = Di), size = 3, alpha = 0.5, colour = "red",position=position_jitter(0.2)) + 
  geom_jitter(data = Y2, aes(x = as.factor(Period), y = Di), size = 3, alpha = 0.5, colour = "cornflowerblue",position=position_jitter(0.2)) + 
  geom_segment(aes(x = 0.5, y = 2.82, xend = 1.5, yend = 2.82))+
  geom_segment(aes(x = 1.5, y = 3.06, xend = 2.5, yend = 3.06))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional distinctiveness")+ scale_x_discrete(labels = c("Historical", "Modern"))

p


#Compute Uniqueness for each species
ui_df = uniqueness_stack(datTOCHE, "Species", dist_mat)

head(ui_df)
write.csv(ui_df, file = "Uniqueness_TOCHE_corr.csv", row.names = FALSE)

### GRaph POP

Ui_Sum <- ui_df %>% 
  group_by(Species) %>%
  summarise_at(vars(Ui), funs(mean(., na.rm=TRUE)))

Ui_Sum <- Ui_Sum[order(Ui_Sum$Ui),] %>% 
  mutate(ID = c(1:239)) %>% 
  left_join(trait, by = "Species")

Ui_Sum <- Ui_Sum %>% 
  filter(!is.na(POP_TOCHE))
########
#Randomizartion to compare means between extirpated and new

#Observed means
A <- subset(Ui_Sum, Ui_Sum$POP_TOCHE == "A")
B <- subset(Ui_Sum, Ui_Sum$POP_TOCHE == "B")
C <- subset(Ui_Sum, Ui_Sum$POP_TOCHE == "C")
meanA <- mean(A$Ui)
meanA #0.615

meanB <- mean(B$Ui)
meanB #0.536

meanC <- mean(C$Ui)
meanC #0.482

obsDifAB <- meanA - meanB #0.078
obsDifBC <- meanB - meanC #0.054

iter    <- 999;          # Total iterations (+1 for observed data = 10k)
diff_AB    <- NULL;# To add difference between groups
diff_BC    <- NULL;# To add difference between groups

N_obs <- 239; # Total number of observations


for(i in 1:iter){   
  reshufled <- Ui_Sum
  reshufled$POP_TOCHE   <- sample(reshufled$POP_TOCHE, N_obs, replace = FALSE);
  
  mean_A <- mean(reshufled %>% filter(POP_TOCHE == "A") %>% pull(Ui))
  mean_B <- mean(reshufled %>% filter(POP_TOCHE == "B") %>% pull(Ui))
  mean_C <- mean(reshufled %>% filter(POP_TOCHE == "C") %>% pull(Ui))
  
  diff_simAB <- mean_A - mean_B
  diff_simBC <- mean_B - mean_C
  
  diff_AB[i] <- diff_simAB
  diff_BC[i] <- diff_simBC
}

#P value Extirpated  -  New
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_AB), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifAB, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_AB >= obsDifAB) + 1;
total_generated   <- length(diff_AB) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.24

#P value New - Shared
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_BC), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifBC, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_BC >= obsDifBC) + 1;
total_generated   <- length(diff_BC) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.26


#Boxplot  Ui POP

windows()
p <- ggplot(Ui_Sum, aes(x=POP_TOCHE, y=Ui, fill = POP_TOCHE)) + 
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(shape=16,position=position_jitter(0.2), colour = "black") +
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), labels = c("Extirpated","New additions","Shared species"))+
  ylim(0,4)+
  geom_text(x = 1.5, y = 3.8, label = "p = 0.24")+
  geom_text(x = 2.5, y = 3.8, label = "p = 0.26")+
  
  geom_segment(aes(x = 1, y = 3.5, xend = 1.9, yend = 3.5))+
  geom_segment(aes(x = 2, y = 3.5, xend = 2.9, yend = 3.5))+
  scale_x_discrete(labels = c("Extirpated","New additions","Shared species"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional Uniqueness")

p


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
tr_mi_z2 <- tr_mi_z2[order(tr_mi_z2$Species),]

tr_mi_z2 <- tr_mi_z2[!duplicated(tr_mi_z2[,1]), ]

tr_mi_z2 <- data.frame(tr_mi_z2[,-1], row.names=tr_mi_z2[,1])


matFUSA <- datFUSA %>% 
  pivot_wider(names_from =  Species, values_from = Abundance) %>% remove_rownames %>%
  column_to_rownames(var = "Period") 

matFUSA[is.na(matFUSA)] <- 0

##############
#FUSAGASUGA #####

matFUSA <- as.matrix(matFUSA)


#Compute distance matrix
dist_mat <- compute_dist_matrix(tr_mi_z_FUSA, metric = "euclidean")

# Compute distinctiveness for each species on each site
di_df = distinctiveness_stack(com_df = datFUSA,  # The site x species table
                              sp_col = "Species",  # Name of the species column
                              com = "Period",  # Name of the community column
                              #abund = "Abundance",  # Relative abundances column (facultative)
                              dist_matrix = dist_mat)  # Functional Distance matrix

head(di_df)
write.csv(di_df, file = "Uniqueness_FUSA_corr.csv", row.names = FALSE)

Di_Sum <- di_df %>% 
  group_by(Species) %>%
  summarise_at(vars(Di), funs(mean(., na.rm=TRUE)))

Di_Sum <- Di_Sum[order(Di_Sum$Di),] %>% 
  mutate(ID = c(1:219)) %>% 
  left_join(trait, by = "Species")

Di_Sum <- Di_Sum %>% 
  filter(!is.na(POP_FUSA))
########
#Randomizartion to compare means between extirpated and new

#Observed means
A <- subset(Di_Sum, Di_Sum$POP_FUSA == "A")
B <- subset(Di_Sum, Di_Sum$POP_FUSA == "B")
C <- subset(Di_Sum, Di_Sum$POP_FUSA == "C")
meanA <- mean(A$Di)
meanA #2.52

meanB <- mean(B$Di)
meanB #2.35

meanC <- mean(C$Di)
meanC #2.10

obsDifAB <- meanA - meanB #0.162
obsDifBC <- meanB - meanC #0.251


iter    <- 999;          # Total iterations (+1 for observed data = 10k)
diff_AB    <- NULL;# To add difference between groups
diff_BC    <- NULL;# To add difference between groups

N_obs <- 219; # Total number of observations


for(i in 1:iter){   
  reshufled <- Di_Sum
  reshufled$POP_FUSA   <- sample(reshufled$POP_FUSA, N_obs, replace = FALSE);
  
  mean_A <- mean(reshufled %>% filter(POP_FUSA == "A") %>% pull(Di))
  mean_B <- mean(reshufled %>% filter(POP_FUSA == "B") %>% pull(Di))
  mean_C <- mean(reshufled %>% filter(POP_FUSA == "C") %>% pull(Di))
  
  diff_simAB <- mean_A - mean_B
  diff_simBC <- mean_B - mean_C
  
  diff_AB[i] <- diff_simAB
  diff_BC[i] <- diff_simBC
}

#P value Extirpated  -  New
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_AB), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifAB, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_AB >= obsDifAB) + 1;
total_generated   <- length(diff_AB) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.23

#P value New - Shared
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_BC), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifBC, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_BC >= obsDifBC) + 1;
total_generated   <- length(diff_BC) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.12


### GRaph POP
windows()
p <- ggplot(Di_Sum, aes(x=ID.x, y=Di, fill = POP_FUSA)) + 
  geom_bar(stat = "identity")+ 
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), breaks=c("A","B","C"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional distinctiveness - Di")

p

#Boxplot  Di POP

windows()
p <- ggplot(Di_Sum, aes(x=POP_FUSA, y=Di, fill = POP_FUSA)) + 
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(shape=16,position=position_jitter(0.2), colour = "black") +
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), labels = c("Extirpated","New additions","Shared species"))+
  scale_x_discrete(labels = c("Extirpated","New additions","Shared species"))+
  #ylim(0,12)+
  #Text for significant differences and p values between group. Replace values with comparison results.
  geom_text(x = 1.5, y = 11, label = "p = 0.23")+
  geom_text(x = 2.5, y = 11, label = "p = 0.12")+
  
  geom_segment(aes(x = 1, y = 10.5, xend = 1.9, yend = 10.5))+
  geom_segment(aes(x = 2, y = 10.5, xend = 2.9, yend = 10.5))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional Distinctiveness")

p

#######GRaph F Distinctiveness by period

Y1 <- subset(di_df, di_df$Period == "HISTORICAL")

Y2 <- subset(di_df, di_df$Period == "MODERN")

meanY1 = mean(Y1$Di)
meanY2 = mean(Y2$Di)

d <- di_df %>% 
  group_by(Period) %>%
  summarise_at(vars(Di), funs(mean(., na.rm=TRUE)))

windows()
p <- ggplot(di_df, aes(x=as.factor(Period), y=Di)) + 
  geom_violin(trim = FALSE, colour = "gray40")+ 
  geom_jitter(shape=16, position=position_jitter(0.2), colour = "gray40") +
  geom_jitter(data = Y1, aes(x = as.factor(Period), y = Di), size = 3, alpha = 0.5, colour = "red",position=position_jitter(0.2)) + 
  geom_jitter(data = Y2, aes(x = as.factor(Period), y = Di), size = 3, alpha = 0.5, colour = "cornflowerblue",position=position_jitter(0.2)) + 
  geom_segment(aes(x = 0.5, y = 2.35, xend = 1.5, yend = 2.35))+
  geom_segment(aes(x = 1.5, y = 2.14, xend = 2.5, yend = 2.14))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional distinctiveness")+ scale_x_discrete(labels = c("Historical", "Modern"))

p


#Compute Uniqueness for each species
ui_df = uniqueness_stack(datFUSA, "Species", dist_mat)

head(ui_df)
write.csv(ui_df, file = "Uniqueness_FUSA_corr.csv", row.names = FALSE)


### GRaph POP

Ui_Sum <- ui_df %>% 
  group_by(Species) %>%
  summarise_at(vars(Ui), funs(mean(., na.rm=TRUE)))

Ui_Sum <- Ui_Sum[order(Ui_Sum$Ui),] %>% 
  mutate(ID = c(1:219)) %>% 
  left_join(trait, by = "Species")


Ui_Sum <- Ui_Sum %>% 
  filter(!is.na(POP_FUSA))
########
#Randomizartion to compare means between extirpated and new

#Observed means
A <- subset(Ui_Sum, Ui_Sum$POP_FUSA == "A")
B <- subset(Ui_Sum, Ui_Sum$POP_FUSA == "B")
C <- subset(Ui_Sum, Ui_Sum$POP_FUSA == "C")
meanA <- mean(A$Ui)
meanA #0.61

meanB <- mean(B$Ui)
meanB #0.49

meanC <- mean(C$Ui)
meanC #0.34

obsDifAB <- meanA - meanB #0.11
obsDifBC <- meanB - meanC #0.15

iter    <- 999;          # Total iterations (+1 for observed data = 10k)
diff_AB    <- NULL;# To add difference between groups
diff_BC    <- NULL;# To add difference between groups

N_obs <- 219; # Total number of observations


for(i in 1:iter){   
  reshufled <- Ui_Sum
  reshufled$POP_FUSA   <- sample(reshufled$POP_FUSA, N_obs, replace = FALSE);
  
  mean_A <- mean(reshufled %>% filter(POP_FUSA == "A") %>% pull(Ui))
  mean_B <- mean(reshufled %>% filter(POP_FUSA == "B") %>% pull(Ui))
  mean_C <- mean(reshufled %>% filter(POP_FUSA == "C") %>% pull(Ui))
  
  diff_simAB <- mean_A - mean_B
  diff_simBC <- mean_B - mean_C
  
  diff_AB[i] <- diff_simAB
  diff_BC[i] <- diff_simBC
}

#P value Extirpated  -  New
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_AB), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifAB, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_AB >= obsDifAB) + 1;
total_generated   <- length(diff_AB) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.16

#P value New - Shared
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_BC), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifBC, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_BC >= obsDifBC) + 1;
total_generated   <- length(diff_BC) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.11

#Boxplot  Ui POP

windows()
p <- ggplot(Ui_Sum, aes(x=POP_FUSA, y=Ui, fill = POP_FUSA)) + 
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(shape=16,position=position_jitter(0.2), colour = "black") +
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), labels = c("Extirpated","New additions","Shared species"))+
  ylim(0,5)+
  geom_text(x = 1.5, y = 4.3, label = "p = 0.16")+
  geom_text(x = 2.5, y = 4.3, label = "p = 0.11")+
  
  geom_segment(aes(x = 1, y = 4, xend = 1.9, yend = 4))+
  geom_segment(aes(x = 2, y = 4, xend = 2.9, yend = 4))+
  scale_x_discrete(labels = c("Extirpated","New additions","Shared species"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional Uniqueness")

p


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
tr_mi_z2 <- tr_mi_z2[order(tr_mi_z2$Species),]

tr_mi_z2 <- tr_mi_z2[!duplicated(tr_mi_z2[,1]), ]

tr_mi_z2 <- data.frame(tr_mi_z2[,-1], row.names=tr_mi_z2[,1])


matHON <- datHON %>% 
  pivot_wider(names_from =  Species, values_from = Abundance) %>% remove_rownames %>%
  column_to_rownames(var = "Period") 

matHON[is.na(matHON)] <- 0

##############
#HONDA

matHON <- as.matrix(matHON)


#Compute distance matrix
dist_mat <- compute_dist_matrix(tr_mi_z_HON, metric = "euclidean")

# Compute distinctiveness for each species on each site
di_df = distinctiveness_stack(com_df = datHON,  # The site x species table
                              sp_col = "Species",  # Name of the species column
                              com = "Period",  # Name of the community column
                              #abund = "Abundance",  # Relative abundances column (facultative)
                              dist_matrix = dist_mat)  # Functional Distance matrix

head(di_df)
write.csv(di_df, file = "Distinctiveness_HON_corr.csv", row.names = FALSE)

Di_Sum <- di_df %>% 
  group_by(Species) %>%
  summarise_at(vars(Di), funs(mean(., na.rm=TRUE)))

Di_Sum <- Di_Sum[order(Di_Sum$Di),] %>% 
  mutate(ID = c(1:229)) %>% 
  left_join(trait, by = "Species")

Di_Sum <- Di_Sum %>% 
  filter(!is.na(POP_HON))
########
#Randomizartion to compare means between extirpated and new

#Observed means
A <- subset(Di_Sum, Di_Sum$POP_HON == "A")
B <- subset(Di_Sum, Di_Sum$POP_HON == "B")
C <- subset(Di_Sum, Di_Sum$POP_HON == "C")
meanA <- mean(A$Di)
meanA #2.63

meanB <- mean(B$Di)
meanB #3.05

meanC <- mean(C$Di)
meanC #2.37

obsDifAB <- meanA - meanB #-0.41
obsDifBC <- meanB - meanC #0.67


iter    <- 999;          # Total iterations (+1 for observed data = 10k)
diff_AB    <- NULL;# To add difference between groups
diff_BC    <- NULL;# To add difference between groups

N_obs <- 229; # Total number of observations


for(i in 1:iter){   
  reshufled <- Di_Sum
  reshufled$POP_HON   <- sample(reshufled$POP_HON, N_obs, replace = FALSE);
  
  mean_A <- mean(reshufled %>% filter(POP_HON == "A") %>% pull(Di))
  mean_B <- mean(reshufled %>% filter(POP_HON == "B") %>% pull(Di))
  mean_C <- mean(reshufled %>% filter(POP_HON == "C") %>% pull(Di))
  
  diff_simAB <- mean_A - mean_B
  diff_simBC <- mean_B - mean_C
  
  diff_AB[i] <- diff_simAB
  diff_BC[i] <- diff_simBC
}

#P value Extirpated  -  New
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_AB), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifAB, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_AB <= obsDifAB) + 1;
total_generated   <- length(diff_AB) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.07

#P value New - Shared
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_BC), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifBC, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_BC >= obsDifBC) + 1;
total_generated   <- length(diff_BC) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.001


### GRaph POP
windows()
p <- ggplot(Di_Sum, aes(x=ID.x, y=Di, fill = POP_HON)) + 
  geom_bar(stat = "identity")+ 
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), breaks=c("A","B","C"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional distinctiveness - Di")

p
#Boxplot  Di POP

windows()
p <- ggplot(Di_Sum, aes(x=POP_HON, y=Di, fill = POP_HON)) + 
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(shape=16,position=position_jitter(0.2), colour = "black") +
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), labels = c("Extirpated","New additions","Shared species"))+
  scale_x_discrete(labels = c("Extirpated","New additions","Shared species"))+
  ylim(0,12)+
  #Text for significant differences and p values between group. Replace values with comparison results.
  geom_text(x = 1.5, y = 11, label = "p = 0.07")+
  geom_text(x = 2.5, y = 11, label = "p = 0.001")+
  
  geom_segment(aes(x = 1, y = 10.5, xend = 1.9, yend = 10.5))+
  geom_segment(aes(x = 2, y = 10.5, xend = 2.9, yend = 10.5))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional Distinctiveness")

p

#######GRaph F Distinctiveness by period

Y1 <- subset(di_df, di_df$Period == "HISTORICAL")

Y2 <- subset(di_df, di_df$Period == "MODERN")

meanY1 = mean(Y1$Di)
meanY2 = mean(Y2$Di)


d <- di_df %>% 
  group_by(Period) %>%
  summarise_at(vars(Di), funs(mean(., na.rm=TRUE)))

windows()
p <- ggplot(di_df, aes(x=as.factor(Period), y=Di)) + 
  geom_violin(trim = FALSE, colour = "gray40")+ 
  geom_jitter(shape=16, position=position_jitter(0.2), colour = "gray40") +
  geom_jitter(data = Y1, aes(x = as.factor(Period), y = Di), size = 3, alpha = 0.5, colour = "red",position=position_jitter(0.2)) + 
  geom_jitter(data = Y2, aes(x = as.factor(Period), y = Di), size = 3, alpha = 0.5, colour = "cornflowerblue",position=position_jitter(0.2)) + 
  geom_segment(aes(x = 0.5, y = 2.35, xend = 1.5, yend = 2.35))+
  geom_segment(aes(x = 1.5, y = 2.82, xend = 2.5, yend = 2.82))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional distinctiveness")+ scale_x_discrete(labels = c("Historical", "Modern"))

p


#Compute Uniqueness for each species
ui_df = uniqueness_stack(datHON, "Species", dist_mat)

head(ui_df)
write.csv(ui_df, file = "Uniqueness_HON_corr.csv", row.names = FALSE)

### GRaph POP

Ui_Sum <- ui_df %>% 
  group_by(Species) %>%
  summarise_at(vars(Ui), funs(mean(., na.rm=TRUE)))

Ui_Sum <- Ui_Sum[order(Ui_Sum$Ui),] %>% 
  mutate(ID = c(1:229)) %>% 
  left_join(trait, by = "Species")

Ui_Sum <- Ui_Sum %>% 
  filter(!is.na(POP_HON))
########
#Randomizartion to compare means between extirpated and new

#Observed means
A <- subset(Ui_Sum, Ui_Sum$POP_HON == "A")
B <- subset(Ui_Sum, Ui_Sum$POP_HON == "B")
C <- subset(Ui_Sum, Ui_Sum$POP_HON == "C")
meanA <- mean(A$Ui)
meanA #0.48

meanB <- mean(B$Ui)
meanB #0.63

meanC <- mean(C$Ui)
meanC #0.45

obsDifAB <- meanA - meanB #-0.15
obsDifBC <- meanB - meanC #0.17

iter    <- 999;          # Total iterations (+1 for observed data = 10k)
diff_AB    <- NULL;# To add difference between groups
diff_BC    <- NULL;# To add difference between groups

N_obs <- 229; # Total number of observations


for(i in 1:iter){   
  reshufled <- Ui_Sum
  reshufled$POP_HON   <- sample(reshufled$POP_HON, N_obs, replace = FALSE);
  
  mean_A <- mean(reshufled %>% filter(POP_HON == "A") %>% pull(Ui))
  mean_B <- mean(reshufled %>% filter(POP_HON == "B") %>% pull(Ui))
  mean_C <- mean(reshufled %>% filter(POP_HON == "C") %>% pull(Ui))
  
  diff_simAB <- mean_A - mean_B
  diff_simBC <- mean_B - mean_C
  
  diff_AB[i] <- diff_simAB
  diff_BC[i] <- diff_simBC
}

#P value Extirpated  -  New
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_AB), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifAB, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_AB <= obsDifAB) + 1;
total_generated   <- length(diff_AB) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.09

#P value New - Shared
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_BC), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifBC, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_BC >= obsDifBC) + 1;
total_generated   <- length(diff_BC) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.006

#Boxplot  Ui POP

windows()
p <- ggplot(Ui_Sum, aes(x=POP_HON, y=Ui, fill = POP_HON)) + 
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(shape=16,position=position_jitter(0.2), colour = "black") +
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), labels = c("Extirpated","New additions","Shared species"))+
  ylim(0,5)+
  geom_text(x = 1.5, y = 4.3, label = "p = 0.09")+
  geom_text(x = 2.5, y = 4.3, label = "p = 0.006")+
  
  geom_segment(aes(x = 1, y = 4, xend = 1.9, yend = 4))+
  geom_segment(aes(x = 2, y = 4, xend = 2.9, yend = 4))+
  scale_x_discrete(labels = c("Extirpated","New additions","Shared species"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional Uniqueness")

p


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
tr_mi_z2 <- tr_mi_z2[order(tr_mi_z2$Species),]

tr_mi_z2 <- tr_mi_z2[!duplicated(tr_mi_z2[,1]), ]

tr_mi_z2 <- data.frame(tr_mi_z2[,-1], row.names=tr_mi_z2[,1])


matSAGU <- datSAGU %>% 
  pivot_wider(names_from =  Species, values_from = Abundance) %>% remove_rownames %>%
  column_to_rownames(var = "Period") 

matSAGU[is.na(matSAGU)] <- 0

##############
#SAN AGUSTIN

matSAGU <- as.matrix(matSAGU)


#Compute distance matrix
dist_mat <- compute_dist_matrix(tr_mi_z_SAGU, metric = "euclidean")

# Compute distinctiveness for each species on each site
di_df = distinctiveness_stack(com_df = datSAGU,  # The site x species table
                              sp_col = "Species",  # Name of the species column
                              com = "Period",  # Name of the community column
                              #abund = "Abundance",  # Relative abundances column (facultative)
                              dist_matrix = dist_mat)  # Functional Distance matrix

head(di_df)
write.csv(di_df, file = "Distinctiveness_SAGU_corr.csv", row.names = FALSE)

Di_Sum <- di_df %>% 
  group_by(Species) %>%
  summarise_at(vars(Di), funs(mean(., na.rm=TRUE)))

Di_Sum <- Di_Sum[order(Di_Sum$Di),] %>% 
  mutate(ID = c(1:254)) %>% 
  left_join(trait, by = "Species")

Di_Sum <- Di_Sum %>% 
  filter(!is.na(POP_SAGU))
########
#Randomizartion to compare means between extirpated and new

#Observed means
A <- subset(Di_Sum, Di_Sum$POP_SAGU == "A")
B <- subset(Di_Sum, Di_Sum$POP_SAGU == "B")
C <- subset(Di_Sum, Di_Sum$POP_SAGU == "C")
meanA <- mean(A$Di)
meanA #2.84

meanB <- mean(B$Di)
meanB #2.78

meanC <- mean(C$Di)
meanC #2.49

obsDifAB <- meanA - meanB #0.06
obsDifBC <- meanB - meanC #0.28


iter    <- 999;          # Total iterations (+1 for observed data = 10k)
diff_AB    <- NULL;# To add difference between groups
diff_BC    <- NULL;# To add difference between groups

N_obs <- 254; # Total number of observations


for(i in 1:iter){   
  reshufled <- Di_Sum
  reshufled$POP_SAGU   <- sample(reshufled$POP_SAGU, N_obs, replace = FALSE);
  
  mean_A <- mean(reshufled %>% filter(POP_SAGU == "A") %>% pull(Di))
  mean_B <- mean(reshufled %>% filter(POP_SAGU == "B") %>% pull(Di))
  mean_C <- mean(reshufled %>% filter(POP_SAGU == "C") %>% pull(Di))
  
  diff_simAB <- mean_A - mean_B
  diff_simBC <- mean_B - mean_C
  
  diff_AB[i] <- diff_simAB
  diff_BC[i] <- diff_simBC
}

#P value Extirpated  -  New
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_AB), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifAB, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_AB >= obsDifAB) + 1;
total_generated   <- length(diff_AB) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.35

#P value New - Shared
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_BC), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifBC, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_BC >= obsDifBC) + 1;
total_generated   <- length(diff_BC) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.04

### GRaph POP
windows()
p <- ggplot(Di_Sum, aes(x=ID.x, y=Di, fill = POP_SAGU)) + 
  geom_bar(stat = "identity")+ 
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), breaks=c("A","B","C"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional distinctiveness - Di")

p

#Boxplot  Di POP

windows()
p <- ggplot(Di_Sum, aes(x=POP_SAGU, y=Di, fill = POP_SAGU)) + 
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(shape=16,position=position_jitter(0.2), colour = "black") +
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), labels = c("Extirpated","New additions","Shared species"))+
  scale_x_discrete(labels = c("Extirpated","New additions","Shared species"))+
  ylim(0,12)+
  #Text for significant differences and p values between group. Replace values with comparison results.
  #Text for significant differences and p values between group. Replace values with comparison results.
  geom_text(x = 1.5, y = 11, label = "p = 0.35")+
  geom_text(x = 2.5, y = 11, label = "p = 0.04")+
  
  geom_segment(aes(x = 1, y = 10.5, xend = 1.9, yend = 10.5))+
  geom_segment(aes(x = 2, y = 10.5, xend = 2.9, yend = 10.5))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional Distinctiveness")

p

#######GRaph F Distinctiveness by period

Y1 <- subset(di_df, di_df$Period == "HISTORICAL")

Y2 <- subset(di_df, di_df$Period == "MODERN")

meanY1 = mean(Y1$Di)
meanY2 = mean(Y2$Di)


d <- di_df %>% 
  group_by(Period) %>%
  summarise_at(vars(Di), funs(mean(., na.rm=TRUE)))

windows()
p <- ggplot(di_df, aes(x=as.factor(Period), y=Di)) + 
  geom_violin(trim = FALSE, colour = "gray40")+ 
  geom_jitter(shape=16, position=position_jitter(0.2), colour = "gray40") +
  geom_jitter(data = Y1, aes(x = as.factor(Period), y = Di), size = 3, alpha = 0.5, colour = "red",position=position_jitter(0.2)) + 
  geom_jitter(data = Y2, aes(x = as.factor(Period), y = Di), size = 3, alpha = 0.5, colour = "cornflowerblue",position=position_jitter(0.2)) + 
  geom_segment(aes(x = 0.5, y = 2.60, xend = 1.5, yend = 2.60))+
  geom_segment(aes(x = 1.5, y = 2.64, xend = 2.5, yend = 2.64))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional distinctiveness")+ scale_x_discrete(labels = c("Historical", "Modern"))

p


#Compute Uniqueness for each species
ui_df = uniqueness_stack(datSAGU, "Species", dist_mat)

head(ui_df)
write.csv(ui_df, file = "Uniqueness_SAGU_corr.csv", row.names = FALSE)

### GRaph POP

Ui_Sum <- ui_df %>% 
  group_by(Species) %>%
  summarise_at(vars(Ui), funs(mean(., na.rm=TRUE)))

Ui_Sum <- Ui_Sum[order(Ui_Sum$Ui),] %>% 
  mutate(ID = c(1:254)) %>% 
  left_join(trait, by = "Species")

Ui_Sum <- Ui_Sum %>% 
  filter(!is.na(POP_SAGU))
########
#Randomizartion to compare means between extirpated and new

#Observed means
A <- subset(Ui_Sum, Ui_Sum$POP_SAGU == "A")
B <- subset(Ui_Sum, Ui_Sum$POP_SAGU == "B")
C <- subset(Ui_Sum, Ui_Sum$POP_SAGU == "C")
meanA <- mean(A$Ui)
meanA #0.62

meanB <- mean(B$Ui)
meanB #0.55

meanC <- mean(C$Ui)
meanC #0.42

obsDifAB <- meanA - meanB #0.06
obsDifBC <- meanB - meanC #0.13

iter    <- 999;          # Total iterations (+1 for observed data = 10k)
diff_AB    <- NULL;# To add difference between groups
diff_BC    <- NULL;# To add difference between groups

N_obs <- 254; # Total number of observations


for(i in 1:iter){   
  reshufled <- Ui_Sum
  reshufled$POP_SAGU   <- sample(reshufled$POP_SAGU, N_obs, replace = FALSE);
  
  mean_A <- mean(reshufled %>% filter(POP_SAGU == "A") %>% pull(Ui))
  mean_B <- mean(reshufled %>% filter(POP_SAGU == "B") %>% pull(Ui))
  mean_C <- mean(reshufled %>% filter(POP_SAGU == "C") %>% pull(Ui))
  
  diff_simAB <- mean_A - mean_B
  diff_simBC <- mean_B - mean_C
  
  diff_AB[i] <- diff_simAB
  diff_BC[i] <- diff_simBC
}

#P value Extirpated  -  New
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_AB), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifAB, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_AB >= obsDifAB) + 1;
total_generated   <- length(diff_AB) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.26

#P value New - Shared
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_BC), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifBC, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_BC >= obsDifBC) + 1;
total_generated   <- length(diff_BC) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.06

#Boxplot  Ui POP

windows()
p <- ggplot(Ui_Sum, aes(x=POP_SAGU, y=Ui, fill = POP_SAGU)) + 
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(shape=16,position=position_jitter(0.2), colour = "black") +
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), labels = c("Extirpated","New additions","Shared species"))+
  ylim(0,5)+
  geom_text(x = 1.5, y = 4.3, label = "p = 0.26")+
  geom_text(x = 2.5, y = 4.3, label = "p = 0.06")+
  
  geom_segment(aes(x = 1, y = 4, xend = 1.9, yend = 4))+
  geom_segment(aes(x = 2, y = 4, xend = 2.9, yend = 4))+
  scale_x_discrete(labels = c("Extirpated","New additions","Shared species"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional Uniqueness")

p


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
tr_mi_z2 <- tr_mi_z2[order(tr_mi_z2$Species),]

tr_mi_z2 <- tr_mi_z2[!duplicated(tr_mi_z2[,1]), ]

tr_mi_z2 <- data.frame(tr_mi_z2[,-1], row.names=tr_mi_z2[,1])


matFLOR <- datFLOR %>% 
  pivot_wider(names_from =  Species, values_from = Abundance) %>% remove_rownames %>%
  column_to_rownames(var = "Period") 

matFLOR[is.na(matFLOR)] <- 0


##############
#FLORENCIA

matFLOR <- as.matrix(matFLOR)


#Compute distance matrix
dist_mat <- compute_dist_matrix(tr_mi_z_FLOR, metric = "euclidean")

# Compute distinctiveness for each species on each site
di_df = distinctiveness_stack(com_df = datFLOR,  # The site x species table
                              sp_col = "Species",  # Name of the species column
                              com = "Period",  # Name of the community column
                              #abund = "Abundance",  # Relative abundances column (facultative)
                              dist_matrix = dist_mat)  # Functional Distance matrix

head(di_df)
write.csv(di_df, file = "Distinctiveness_FLOR_corr.csv", row.names = FALSE)

Di_Sum <- di_df %>% 
  group_by(Species) %>%
  summarise_at(vars(Di), funs(mean(., na.rm=TRUE)))

Di_Sum <- Di_Sum[order(Di_Sum$Di),] %>% 
  mutate(ID = c(1:287)) %>% 
  left_join(trait, by = "Species")

Di_Sum <- Di_Sum %>% 
  filter(!is.na(POP_FLOR))
########
#Randomizartion to compare means between extirpated and new

#Observed means
A <- subset(Di_Sum, Di_Sum$POP_FLOR == "A")
B <- subset(Di_Sum, Di_Sum$POP_FLOR == "B")
C <- subset(Di_Sum, Di_Sum$POP_FLOR == "C")
meanA <- mean(A$Di)
meanA #3.00

meanB <- mean(B$Di)
meanB #3.82

meanC <- mean(C$Di)
meanC #3.18

obsDifAB <- meanA - meanB #-0.81
obsDifBC <- meanB - meanC #0.63


iter    <- 999;          # Total iterations (+1 for observed data = 10k)
diff_AB    <- NULL;# To add difference between groups
diff_BC    <- NULL;# To add difference between groups

N_obs <- 287; # Total number of observations


for(i in 1:iter){   
  reshufled <- Di_Sum
  reshufled$POP_FLOR   <- sample(reshufled$POP_FLOR, N_obs, replace = FALSE);
  
  mean_A <- mean(reshufled %>% filter(POP_FLOR == "A") %>% pull(Di))
  mean_B <- mean(reshufled %>% filter(POP_FLOR == "B") %>% pull(Di))
  mean_C <- mean(reshufled %>% filter(POP_FLOR == "C") %>% pull(Di))
  
  diff_simAB <- mean_A - mean_B
  diff_simBC <- mean_B - mean_C
  
  diff_AB[i] <- diff_simAB
  diff_BC[i] <- diff_simBC
}

#P value Extirpated  -  New
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_AB), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifAB, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_AB <= obsDifAB) + 1;
total_generated   <- length(diff_AB) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.001

#P value New - Shared
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_BC), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifBC, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_BC >= obsDifBC) + 1;
total_generated   <- length(diff_BC) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.002

### GRaph POP
windows()
p <- ggplot(Di_Sum, aes(x=ID.x, y=Di, fill = POP_FLOR)) + 
  geom_bar(stat = "identity")+ 
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), breaks=c("A","B","C"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional distinctiveness - Di")

p
#Boxplot  Di POP

windows()
p <- ggplot(Di_Sum, aes(x=POP_FLOR, y=Di, fill = POP_FLOR)) + 
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(shape=16,position=position_jitter(0.2), colour = "black") +
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), labels = c("Extirpated","New additions","Shared species"))+
  scale_x_discrete(labels = c("Extirpated","New additions","Shared species"))+
  #ylim(0,12)+
  #Text for significant differences and p values between group. Replace values with comparison results.
  geom_text(x = 1.5, y = 11, label = "p = 0.001")+
  geom_text(x = 2.5, y = 11, label = "p = 0.002")+
  
  geom_segment(aes(x = 1, y = 10.5, xend = 1.9, yend = 10.5))+
  geom_segment(aes(x = 2, y = 10.5, xend = 2.9, yend = 10.5))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional Distinctiveness")

p

#######GRaph F Distinctiveness by period

Y1 <- subset(di_df, di_df$Period == "HISTORICAL")

Y2 <- subset(di_df, di_df$Period == "MODERN")

meanY1 = mean(Y1$Di)
meanY2 = mean(Y2$Di)

d <- di_df %>% 
  group_by(Period) %>%
  summarise_at(vars(Di), funs(mean(., na.rm=TRUE)))

windows()
p <- ggplot(di_df, aes(x=as.factor(Period), y=Di)) + 
  geom_violin(trim = FALSE, colour = "gray40")+ 
  geom_jitter(shape=16, position=position_jitter(0.2), colour = "gray40") +
  geom_jitter(data = Y1, aes(x = as.factor(Period), y = Di), size = 3, alpha = 0.5, colour = "red",position=position_jitter(0.2)) + 
  geom_jitter(data = Y2, aes(x = as.factor(Period), y = Di), size = 3, alpha = 0.5, colour = "cornflowerblue",position=position_jitter(0.2)) + 
  geom_segment(aes(x = 0.5, y = 3.01, xend = 1.5, yend = 3.01))+
  geom_segment(aes(x = 1.5, y = 3.61, xend = 2.5, yend = 3.61))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional distinctiveness")+ scale_x_discrete(labels = c("Historical", "Modern"))

p


#Compute Uniqueness for each species

ui_df = uniqueness_stack(datFLOR, "Species", dist_mat)

head(ui_df)
write.csv(ui_df, file = "Uniqueness_FLOR_corr.csv", row.names = FALSE)

### GRaph POP

Ui_Sum <- ui_df %>% 
  group_by(Species) %>%
  summarise_at(vars(Ui), funs(mean(., na.rm=TRUE)))

Ui_Sum <- Ui_Sum[order(Ui_Sum$Ui),] %>% 
  mutate(ID = c(1:287)) %>% 
  left_join(trait, by = "Species")

Ui_Sum <- Ui_Sum %>% 
  filter(!is.na(POP_FLOR))
########
#Randomizartion to compare means between extirpated and new

#Observed means
A <- subset(Ui_Sum, Ui_Sum$POP_FLOR == "A")
B <- subset(Ui_Sum, Ui_Sum$POP_FLOR == "B")
C <- subset(Ui_Sum, Ui_Sum$POP_FLOR == "C")
meanA <- mean(A$Ui)
meanA #0.58

meanB <- mean(B$Ui)
meanB #0.67

meanC <- mean(C$Ui)
meanC #0.51

obsDifAB <- meanA - meanB #-0.09
obsDifBC <- meanB - meanC #0.15

iter    <- 999;          # Total iterations (+1 for observed data = 10k)
diff_AB    <- NULL;# To add difference between groups
diff_BC    <- NULL;# To add difference between groups

N_obs <- 287; # Total number of observations


for(i in 1:iter){   
  reshufled <- Ui_Sum
  reshufled$POP_FLOR   <- sample(reshufled$POP_FLOR, N_obs, replace = FALSE);
  
  mean_A <- mean(reshufled %>% filter(POP_FLOR == "A") %>% pull(Ui))
  mean_B <- mean(reshufled %>% filter(POP_FLOR == "B") %>% pull(Ui))
  mean_C <- mean(reshufled %>% filter(POP_FLOR == "C") %>% pull(Ui))
  
  diff_simAB <- mean_A - mean_B
  diff_simBC <- mean_B - mean_C
  
  diff_AB[i] <- diff_simAB
  diff_BC[i] <- diff_simBC
}

#P value Extirpated  -  New
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_AB), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifAB, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_AB <= obsDifAB) + 1;
total_generated   <- length(diff_AB) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.13

#P value New - Shared
windows()
ggplot() +
  ylab("Count") + xlab("Simulated mean difference") +
  geom_histogram(aes(x = diff_BC), bins = 30, 
                 fill = "grey", alpha = 0.4, colour = "black") +
  geom_vline(xintercept = obsDifBC, linewidth = 1, 
             linetype = "dashed", colour = "black") + 
  theme_classic()


less_or_equal_obs <- sum(diff_BC >= obsDifBC) + 1;
total_generated   <- length(diff_BC) + 1;
new_p_value       <- less_or_equal_obs / total_generated;
# p = 0.05

#Boxplot  Ui POP

windows()
p <- ggplot(Ui_Sum, aes(x=POP_FLOR, y=Ui, fill = POP_FLOR)) + 
  geom_boxplot(outlier.shape = NA)+ 
  geom_jitter(shape=16,position=position_jitter(0.2), colour = "black") +
  scale_fill_manual(values=c("red", "cornflowerblue", "azure3"), labels = c("Extirpated","New additions","Shared species"))+
  #ylim(0,3)+
  geom_text(x = 1.5, y = 4.3, label = "p = 0.13")+
  geom_text(x = 2.5, y = 4.3, label = "p = 0.05")+
  
  geom_segment(aes(x = 1, y = 4, xend = 1.9, yend = 4))+
  geom_segment(aes(x = 2, y = 4, xend = 2.9, yend = 4))+
  scale_x_discrete(labels = c("Extirpated","New additions","Shared species"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  labs(x = "", y = "Functional Uniqueness")

p


