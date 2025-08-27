#1. Functional space and hypervolume estimation by region
library(pacman)
pacman::p_load(dplyr, plyr, readr, tbible, FD, ade4, cowplot, mice, reshape2, tidyr, ks, hypervolume, alphallhu, purrr, TTR, plotrix, agricolae, psych)

library(factoextra)
library(ggrepel)
library(tibble)
library(reshape2)

#### Load data

# load: trait data
trait <- read.csv("C:\\Users\\camil\\OneDrive - SELVA Investigación para la conservación en el neotropico\\Documentos\\Camila\\Publicaciones\\Chapman_Historical Occupancy\\All_species_traits_CorrJuli_FINAL_spLOST.csv", sep = ";")

head(trait)
#Remove new species

trait <- trait %>% 
  filter(!HISTORICAL %in% c(0, NA))

# Exploring trait correlations (verify trait column numbers)

just_trait <- trait[, c(67,68,69,70,71,72,73,74,75,76,77)]
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
  select(Species,MODERN, HISTORICAL,ANT_H,ANT_M, ANT_EXT,MIR_H_corr,MIR_M_LOST,MIR_EXT,BARB_H_corr,BARB_M_LOST,BARB_EXT,SAGU_H_corr, SAGU_M_LOST,SAGU_EXT,TOCHE_H_corr,TOCHE_M_LOST,TOCHE_EXT, HON_H_corr,HON_M_LOST,HON_EXT, FLOR_H_corr, FLOR_M_LOST,FLOR_EXT, FUSA_H_corr, FUSA_M_LOST,FUSA_EXT, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) 

ANT <- trait_z %>% 
  filter(ANT_H >= 0)

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
## PCA ## SAN ANTONIO
#### --------------------------------------------------------------

#table just with columns of morphological traits and species names as rownames
tr_mi_z_ANT <- ANT %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)%>% 
  column_to_rownames(var = "Species")

#Run PCA
pcaTotal_ANT <- princomp(tr_mi_z_ANT, cor = TRUE, scores = TRUE)
summary(pcaTotal_ANT)

#Scores
scoresPCATotal_ANT <- as.data.frame(pcaTotal_ANT$scores) %>% 
  tibble::rownames_to_column("Species")

scoresPCATotal_ANT <- scoresPCATotal_ANT %>% 
  # convert long to wide
  tidyr::gather(key, value, -Species) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# loadings
loadingsPCATotal_ANT <- as.data.frame(unclass(pcaTotal_ANT$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()

loadingsPCATotal_ANT <- loadingsPCATotal_ANT %>% 
  # convert long to wide
  tidyr::gather(key, value, -trait) %>% 
  tidyr::unite(col, key) %>% 
  tidyr::spread(col, value) 

# scalar to adjust arrow size
sc_mi <- 7

loadings_sc_ANT <- loadingsPCATotal_ANT %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")), funs(.*sc_mi)) %>% 
  # variable names
  mutate(trait = c("Bill width","Hand-Wing Index", "Total culmen", "Tail length", "Tarsus", "Body mass", "Wing length"))

#######################


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

ANT_H <- trait_z %>% 
  filter(ANT_H >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_ANT, by = "Species") 

ANT_M <- trait_z %>% 
  filter(ANT_M >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_ANT, by = "Species")

MIR_H <- trait_z %>% 
  filter(MIR_H_corr >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_MIR, by = "Species") 

MIR_M <- trait_z %>% 
  filter(MIR_M_LOST >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_MIR, by = "Species")

BARB_H <- trait_z %>% 
  filter(BARB_H_corr >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_BARB, by = "Species") 

BARB_M <- trait_z %>% 
  filter(BARB_M_LOST >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_BARB, by = "Species")

SAGU_H <- trait_z %>% 
  filter(SAGU_H_corr >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_SAGU, by = "Species") 

SAGU_M <- trait_z %>% 
  filter(SAGU_M_LOST >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_SAGU, by = "Species")

TOCHE_H <- trait_z %>% 
  filter(TOCHE_H_corr >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_TOCHE, by = "Species") 

TOCHE_M <- trait_z %>% 
  filter(TOCHE_M_LOST >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_TOCHE, by = "Species")

HON_H <- trait_z %>% 
  filter(HON_H_corr >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_HON, by = "Species") 

HON_M <- trait_z %>% 
  filter(HON_M_LOST >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_HON, by = "Species")

FLOR_H <- trait_z %>% 
  filter(FLOR_H_corr >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FLOR, by = "Species") 

FLOR_M <- trait_z %>% 
  filter(FLOR_M_LOST >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FLOR, by = "Species")

FUSA_H <- trait_z %>% 
  filter(FUSA_H_corr >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FUSA, by = "Species") 

FUSA_M <- trait_z %>% 
  filter(FUSA_M_LOST >= 1) %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FUSA, by = "Species")

# kernel density estimation for each Region and each period

pc_raw_ANT_H <- ANT_H %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "Species")

pc_raw_ANT_M <- ANT_M %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "Species")

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
hpi_ANT_H <- Hpi(x = pc_raw_ANT_H)
hpi_ANT_M <- Hpi(x = pc_raw_ANT_M)
 
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

# kernel density estimation San Antonio_Historical  
est_ANT_H <- kde(x = pc_raw_ANT_H, H = hpi_ANT_H, compute.cont = TRUE)  

den_ANT_H <- list(est_ANT_H$eval.points[[1]], est_ANT_H$eval.points[[2]], est_ANT_H$estimate)
names(den_ANT_H) <- c("x", "y", "z")
dimnames(den_ANT_H$z) <- list(den_ANT_H$x, den_ANT_H$y)
dcc_ANT_H <- reshape2::melt(den_ANT_H$z)

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
cl_50_ANT_H <- cl(df = den_ANT_H, prob = 0.50)
# 0.95 probability kernel
cl_95_ANT_H <- cl(df = den_ANT_H, prob = 0.95)
# 0.99 probability kernel
cl_99_ANT_H <- cl(df = den_ANT_H, prob = 0.99)

# kernel density estimation San Antonio Modern

est_ANT_M <- kde(x = pc_raw_ANT_M, H = hpi_ANT_M, compute.cont = TRUE)  

den_ANT_M <- list(est_ANT_M$eval.points[[1]], est_ANT_M$eval.points[[2]], est_ANT_M$estimate)
names(den_ANT_M) <- c("x", "y", "z")
dimnames(den_ANT_M$z) <- list(den_ANT_M$x, den_ANT_M$y)
dcc_ANT_M <- reshape2::melt(den_ANT_M$z)

# 0.5 probability kernel
cl_50_ANT_M <- cl(df = den_ANT_M, prob = 0.50)
# 0.95 probability kernel
cl_95_ANT_M <- cl(df = den_ANT_M, prob = 0.95)
# 0.99 probability kernel
cl_99_ANT_M <- cl(df = den_ANT_M, prob = 0.99)


# kernel density estimation Miraflores_Historical  
est_MIR_H <- kde(x = pc_raw_MIR_H, H = hpi_MIR_H, compute.cont = TRUE)  

den_MIR_H <- list(est_MIR_H$eval.points[[1]], est_MIR_H$eval.points[[2]], est_MIR_H$estimate)
names(den_MIR_H) <- c("x", "y", "z")
dimnames(den_MIR_H$z) <- list(den_MIR_H$x, den_MIR_H$y)
dcc_MIR_H <- reshape2::melt(den_MIR_H$z)

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
dcc_MIR_M <- reshape2::melt(den_MIR_M$z)

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
dcc_BARB_H <- reshape2::melt(den_BARB_H$z)

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
dcc_BARB_M <- reshape2::melt(den_BARB_M$z)

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
dcc_SAGU_M <- reshape2::melt(den_SAGU_M$z)

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
dcc_SAGU_H <- reshape2::melt(den_SAGU_H$z)

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
dcc_TOCHE_M <- reshape2::melt(den_TOCHE_M$z)

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
dcc_TOCHE_H <- reshape2::melt(den_TOCHE_H$z)

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
dcc_HON_M <- reshape2::melt(den_HON_M$z)

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
dcc_HON_H <- reshape2::melt(den_HON_H$z)

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
dcc_FLOR_M <- reshape2::melt(den_FLOR_M$z)

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
dcc_FLOR_H <- reshape2::melt(den_FLOR_H$z)

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
dcc_FUSA_M <- reshape2::melt(den_FUSA_M$z)

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
dcc_FUSA_H <- reshape2::melt(den_FUSA_H$z)

# 0.5 probability kernel
cl_50_FUSA_H <- cl(df = den_FUSA_H, prob = 0.50)
# 0.95 probability kernel
cl_95_FUSA_H <- cl(df = den_FUSA_H, prob = 0.95)
# 0.99 probability kernel
cl_99_FUSA_H <- cl(df = den_FUSA_H, prob = 0.99)

# save principal component data
PCA_ANT <- trait_z %>% 
  filter(ANT_H >= 0) %>% 
  select(Species,ANT_H, ANT_M,ANT_EXT, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_ANT, by = "Species")

PCA_MIR <- trait_z %>% 
  filter(MIR_H_corr >= 0) %>% 
  select(Species,MIR_H_corr, MIR_M_LOST,MIR_EXT,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_MIR, by = "Species") 

PCA_BARB <- trait_z %>% 
  filter(BARB_H_corr >= 0) %>% 
  select(Species,BARB_H_corr, BARB_M_LOST,BARB_EXT, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>%  
  left_join(scoresPCATotal_BARB, by = "Species") 

PCA_SAGU <- trait_z %>% 
  filter(SAGU_H_corr >= 0) %>% 
  select(Species,SAGU_H_corr, SAGU_M_LOST,SAGU_EXT, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_SAGU, by = "Species") 

PCA_TOCHE <- trait_z %>% 
  filter(TOCHE_H_corr >= 0) %>% 
  select(Species,TOCHE_H_corr, TOCHE_M_LOST,TOCHE_EXT, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_TOCHE, by = "Species") 

PCA_HON <- trait_z %>% 
  filter(HON_H_corr >= 0) %>% 
  select(Species,HON_H_corr,HON_M_LOST,HON_EXT, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_HON, by = "Species") 

PCA_FLOR <- trait_z %>% 
  filter(FLOR_H_corr >= 0) %>% 
  select(Species,FLOR_H_corr, FLOR_M_LOST,FLOR_EXT, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FLOR, by = "Species") 

PCA_FUSA <- trait_z %>% 
  filter(FUSA_H_corr >= 0) %>% 
  select(Species,FUSA_H_corr,FUSA_M_LOST,FUSA_EXT, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  left_join(scoresPCATotal_FUSA, by = "Species") 

write.csv(PCA_ANT, file = "PCA_ANT_corr_FINAL_LOST.csv", row.names = FALSE)
write.csv(PCA_MIR, file = "PCA_MIR_corr_FINAL_LOST.csv", row.names = FALSE)
write.csv(PCA_BARB, file = "PCA_BARB_corr_FINAL_LOST.csv", row.names = FALSE)
write.csv(PCA_SAGU, file = "PCA_SAGU_corr_FINAL_LOST.csv", row.names = FALSE)
write.csv(PCA_TOCHE, file = "PCA_TOCHE_corr_FINAL_LOST.csv", row.names = FALSE)
write.csv(PCA_HON, file = "PCA_HON_corr_FINAL_LOST.csv", row.names = FALSE)
write.csv(PCA_FLOR, file = "PCA_FLOR_corr_FINAL_LOST.csv", row.names = FALSE)
write.csv(PCA_FUSA, file = "PCA_FUSA_corr_FINAL_LOST.csv", row.names = FALSE)

write.csv(scoresPCATotal_ANT, file = "PCA_Total_ANT_corr_FINAL_LOST.csv")
write.csv(scoresPCATotal_MIR, file = "PCA_Total_MIR_corr_FINAL_LOST.csv")
write.csv(scoresPCATotal_BARB, file = "PCA_Total_BARB_corr_FINAL_LOST.csv")
write.csv(scoresPCATotal_TOCHE, file = "PCA_Total_TOCHE_corr_FINAL_LOST.csv")
write.csv(scoresPCATotal_SAGU, file = "PCA_Total_SAGU_corr_FINAL_LOST.csv")
write.csv(scoresPCATotal_HON, file = "PCA_Total_HON_corr_FINAL_LOST.csv")
write.csv(scoresPCATotal_FLOR, file = "PCA_Total_FLOR_corr_FINAL_LOST.csv")
write.csv(scoresPCATotal_FUSA, file = "PCA_Total_FUSA_corr_FINAL_LOST.csv")

write.csv(loadingsPCATotal_ANT, file = "Loadings_Total_ANT_corr_FINAL_LOST.csv", row.names = FALSE)
write.csv(loadingsPCATotal_MIR, file = "Loadings_Total_MIR_corr_FINAL_LOST.csv", row.names = FALSE)
write.csv(loadingsPCATotal_BARB, file = "Loadings_Total_BARB_corr_FINAL_LOST.csv", row.names = FALSE)
write.csv(loadingsPCATotal_SAGU, file = "Loadings_Total_SAGU.csv_corr_FINAL_LOST", row.names = FALSE)
write.csv(loadingsPCATotal_TOCHE, file = "Loadings_Total_TOCHE_corr_FINAL_LOST.csv", row.names = FALSE)
write.csv(loadingsPCATotal_HON, file = "Loadings_Total_HON_corr_FINAL_LOST.csv", row.names = FALSE)
write.csv(loadingsPCATotal_FLOR, file = "Loadings_Total_FLOR_corr_FINAL_LOST.csv", row.names = FALSE)
write.csv(loadingsPCATotal_FUSA, file = "Loadings_Total_FUSA_corr_FINAL_LOST.csv", row.names = FALSE)

#### PCA plots
library(ggplot2)

# colour palette
col_pal <- colorRampPalette(c("grey30","grey60", "grey80", "white"))(200)


### Folowing lines indentify species that correspond to extinction,
#recolonization or new additions in different periods

extir_ANT <- subset(PCA_ANT, ANT_EXT == 1)

extir_MIR <- subset(PCA_MIR, MIR_EXT == 1)
pextir_MIR <- subset(PCA_MIR, MIR_EXT == 2)

extir_BARB <- subset(PCA_BARB, BARB_EXT == 1)
pextir_BARB <- subset(PCA_BARB, BARB_EXT == 2)

extir_SAGU <- subset(PCA_SAGU, SAGU_EXT == 1)
pextir_SAGU <- subset(PCA_MIR, MIR_EXT == 2)

extir_TOCHE <- subset(PCA_TOCHE, TOCHE_EXT == 1)
pextir_TOCHE <- subset(PCA_TOCHE, TOCHE_EXT == 2)

extir_HON <- subset(PCA_HON, HON_EXT == 1)
pextir_HON <- subset(PCA_HON, HON_EXT == 2)

extir_FLOR <- subset(PCA_FLOR, FLOR_EXT == 1)
pextir_FLOR <- subset(PCA_FLOR, FLOR_EXT == 2)

extir_FUSA <- subset(PCA_FUSA, FUSA_EXT == 1)
pextir_FUSA <- subset(PCA_FUSA, FUSA_EXT == 2)

# plot San Antonio Historical
pca_plot_ANT_H <- ggplot(dcc_ANT_H, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_ANT, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_ANT, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  
  # add arrows
  geom_segment(data = loadings_sc_ANT, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_ANT_H, colour = "grey30", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_95_ANT_H, colour = "grey60", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_99_ANT_H, colour = "grey70", linewidth = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_ANT, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_ANT, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_ANT, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
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
  geom_text(x= 5, y=8, label="San Antonio Historical", color="black", size = 5)

# display plot
windows()
pca_plot_ANT_H

# plot SAN ANTONIO Modern
pca_plot_ANT_M <- ggplot(dcc_ANT_M, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_ANT, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_ANT, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  
  # add arrows
  geom_segment(data = loadings_sc_ANT, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_ANT_M, colour = "grey30", size = 1) +
  geom_contour(aes(z = value), breaks = cl_95_ANT_M, colour = "grey60", size = 1) +
  geom_contour(aes(z = value), breaks = cl_99_ANT_M, colour = "grey70", size = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_ANT, aes(x = 0, y = 0, xend = Comp.1, yend = Comp.2), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_ANT, aes(x = 0, y = 0, xend = -Comp.1, yend = -Comp.2), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_ANT, aes(x = Comp.1, y = Comp.2, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
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
  geom_text(x= 5, y=8, label="San Antonio Modern", color="black", size = 5)

# display plot
windows()
pca_plot_ANT_M


# plot Miraflores Historical
pca_plot_MIR_H <- ggplot(dcc_MIR_H, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_MIR, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_MIR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = pextir_MIR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  labs(x = "PC1 - Size (66%)", y = "PC2 - Dispersal ability (21%)") +
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
  geom_point(data = pextir_MIR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
 
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
  labs(x = "PC1 - Size (66%)", y = "PC2 - Dispersal ability (21%)") +
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
  geom_point(data = pextir_BARB, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  labs(x = "PC1 - Size (59%)", y = "PC2 - Dispersal ability (22%)") +
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
  geom_point(data = pextir_BARB, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  labs(x = "PC1 - Size (59%)", y = "PC2 - Dispersal ability (22%)") +
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
  geom_point(data = pextir_SAGU, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  labs(x = "PC1 - Size (66%)", y = "PC2 - Dispersal ability (21%)") +
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
  geom_point(data = pextir_SAGU, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  labs(x = "PC1 - Size (66%)", y = "PC2 - Dispersal ability (21%)") +
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
  geom_point(data = pextir_TOCHE, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  labs(x = "PC1 - Size (65%)", y = "PC2 - Dispersal ability (21%)") +
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
  geom_point(data = pextir_TOCHE, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  labs(x = "PC1 - Size (65%)", y = "PC2 - Dispersal ability (21%)") +
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
  geom_point(data = pextir_HON, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  geom_point(data = pextir_HON, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  geom_point(data = pextir_FLOR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  geom_point(data = pextir_FLOR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  geom_point(data = pextir_FUSA, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  geom_point(data = pextir_FUSA, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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

#SAN ANTONIO
PCA_ANT_H_A <- subset(PCA_ANT, ANT_H >= 1)

PCA_ANT_H_A <-  PCA_ANT_H_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

PCA_ANT_M_A <- subset(PCA_ANT, ANT_M >= 1)

PCA_ANT_M_A <-  PCA_ANT_M_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")


#MIRAFLORES
PCA_MIR_H_A <- subset(PCA_MIR, MIR_H_corr >= 1)

PCA_MIR_H_A <-  PCA_MIR_H_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

PCA_MIR_M_A <- subset(PCA_MIR, MIR_M_LOST >= 1)

PCA_MIR_M_A <-  PCA_MIR_M_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

#BARBACOAS
PCA_BARB_H_A <- subset(PCA_BARB, BARB_H_corr >= 1)

PCA_BARB_H_A <-  PCA_BARB_H_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

PCA_BARB_M_A <- subset(PCA_BARB, BARB_M_LOST >= 1)

PCA_BARB_M_A <-  PCA_BARB_M_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

#TOCHE
PCA_TOCHE_H_A <- subset(PCA_TOCHE, TOCHE_H_corr >= 1)

PCA_TOCHE_H_A <-  PCA_TOCHE_H_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species") 

PCA_TOCHE_M_A <- subset(PCA_TOCHE, TOCHE_M_LOST >= 1)

PCA_TOCHE_M_A <-  PCA_TOCHE_M_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

#SAN AGUSTIN
PCA_SAGU_H_A <- subset(PCA_SAGU, SAGU_H_corr >= 1)

PCA_SAGU_H_A <-  PCA_SAGU_H_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species") 

PCA_SAGU_M_A <- subset(PCA_SAGU, SAGU_M_LOST >= 1)

PCA_SAGU_M_A <-  PCA_SAGU_M_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

#HONDA
PCA_HON_H_A <- subset(PCA_HON, HON_H_corr >= 1)

PCA_HON_H_A <-  PCA_HON_H_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species") 

PCA_HON_M_A <- subset(PCA_HON, HON_M_LOST >= 1)

PCA_HON_M_A <-  PCA_HON_M_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

#FLORENCIA
PCA_FLOR_H_A <- subset(PCA_FLOR, FLOR_H_corr >= 1)

PCA_FLOR_H_A <-  PCA_FLOR_H_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species") 

PCA_FLOR_M_A <- subset(PCA_FLOR, FLOR_M_LOST >= 1)

PCA_FLOR_M_A <-  PCA_FLOR_M_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

#FUSAGASUGA
PCA_FUSA_H_A <- subset(PCA_FUSA, FUSA_H_corr >= 1)

PCA_FUSA_H_A <-  PCA_FUSA_H_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3)  %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")

PCA_FUSA_M_A <- subset(PCA_FUSA, FUSA_M_LOST >= 1)

PCA_FUSA_M_A <-  PCA_FUSA_M_A %>% 
  select(Species, Comp.1, Comp.2, Comp.3) %>% 
  remove_rownames %>% 
  column_to_rownames(var = "Species")


# hypervolumes: Takes a long time to run if all components included. To test, first run using just 3 components by changing code above to:
# PCA_BARB_H_A <-  PCA_BARB_H_A %>% 
# select(Species, Comp.1, Comp.2, Comp.3) %>% 
# remove_rownames %>% 
# column_to_rownames(var = "Species")


library(hypervolume)

#SAN ANTONIO

set.seed(3)
ANT_H_hyper_mi <- hypervolume::hypervolume_svm(PCA_ANT_H_A, name = "ANT_H")

saveRDS(ANT_H_hyper_mi, "ANT_H_hyper_mi5_corr_FINAL_LOST.rds")

set.seed(3)
ANT_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_ANT_M_A, name = "ANT_M")

saveRDS(ANT_M_hyper_mi, "ANT_M_hyper_mi5_corr_FINAL_LOST.rds")

# load: hyper volumes
#ANT_H_hyper_mi <- readRDS("ANT_H_hyper_mi5.rds")
#ANT_M_hyper_mi <- readRDS("ANT_M_hyper_mi5.rds")


# set hypervolume comparisons 

set.seed(3)
bm_set_mi_ANTHvsANTM <- hypervolume::hypervolume_set(ANT_H_hyper_mi, ANT_M_hyper_mi, check.memory = FALSE)

# overlap statistics
bm_over_mi <- hypervolume::hypervolume_overlap_statistics(bm_set_mi_ANTHvsANTM)

# summarise volumes Historical  vs Modern

hypervolume::get_volume(bm_set_mi_ANTHvsANTM)

# Plot 3D hypervolumes 
library(rgl)
library(alphahull)

#Historical vs Modern
bm_set_mi_ANTHvsANTM@HVList$HV1 <- NULL
bm_set_mi_ANTHvsANTM@HVList$HV2 <- NULL
bm_set_mi_ANTHvsANTM@HVList$Union <- NULL

# trait names
names <- c("PC1", "PC2","PC3")

colnames(bm_set_mi_ANTHvsANTM@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_ANTHvsANTM@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_ANTHvsANTM@HVList$Unique_2@RandomPoints) <- names

# plot hypervolumes 
hypervolume::plot.HypervolumeList(bm_set_mi_ANTHvsANTM, show.3d = TRUE, plot.3d.axes.id = c(1,2,3), show.random = FALSE, show.data = TRUE, show.centroid = FALSE, show.density = TRUE, cex.random = 5, colors = c("azure3", "red","cornflowerblue"))

#Fancy plot 3D Hypervolume
######################################
#Data frame for Plotly Hyper Volume Plot
inter <- data.frame(bm_set_mi_ANTHvsANTM@HVList$Intersection@RandomPoints)
inter <- inter %>% 
  mutate(section = "intersection")

uni1 <- data.frame(bm_set_mi_ANTHvsANTM@HVList$Unique_1@RandomPoints) 
uni1 <- uni1 %>% 
  mutate(section = "Unique Historical")

uni2 <- data.frame(bm_set_mi_ANTHvsANTM@HVList$Unique_2@RandomPoints)
uni2 <- uni2 %>% 
  mutate(section = "Unique Modern")

HV <- rbind(inter, uni1, uni2)
write.csv(HV,"HyperVolume_ANT.csv")

colnames(bm_set_mi_ANTHvsANTM@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_ANTHvsANTM@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_ANTHvsANTM@HVList$Unique_2@RandomPoints) <- names


library(plotly)

HV <- read.csv("HyperVolume_ANT.csv", h = T)


pal <- c("azure3","red", "cornflowerblue" )


fig <- plot_ly(HV, x = ~PC1, y = ~(PC2), z = ~(PC3), type = 'scatter3d', mode = 'markers', color = ~section, colors = pal, marker = list(size = 3,opacity = 0.6, line = list(color = 'grey66', width = 2)))

fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                   yaxis = list(title = 'PC2'),
                                   zaxis = list(title = 'PC3')))

fig

######################################
#Miraflores
set.seed(3)
MIR_H_hyper_mi <- hypervolume::hypervolume_svm(PCA_MIR_H_A, name = "MIR_H")

saveRDS(MIR_H_hyper_mi, "MIR_H_hyper_mi5_corr_FINAL_LOST.rds")

set.seed(3)
MIR_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_MIR_M_A, name = "MIR_M")

saveRDS(MIR_M_hyper_mi, "MIR_M_hyper_mi5_corr_FINAL_LOST.rds")

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
names <- c("PC1", "PC2","PC3")

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


fig <- plot_ly(HV, x = ~PC1, y = ~(PC2), z = ~(PC3), type = 'scatter3d', mode = 'markers', color = ~section, colors = pal, marker = list(size = 3,opacity = 0.6, line = list(color = 'grey66', width = 2)))

fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                   yaxis = list(title = 'PC2'),
                                   zaxis = list(title = 'PC3')))

fig

######################################

#Barbacoas 
set.seed(3)
BARB_H_hyper_mi <- hypervolume::hypervolume_svm(PCA_BARB_H_A, name = "BARB_H")

saveRDS(BARB_H_hyper_mi, "BARB_H_hyper_mi5_corr_FINAL_LOST.rds")

set.seed(3)
BARB_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_BARB_M_A, name = "BARB_M")

saveRDS(BARB_M_hyper_mi, "BARB_M_hyper_mi5_corr_FINAL_LOST.rds")

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
names <- c("PC1", "PC2","PC3")

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

saveRDS(TOCHE_H_hyper_mi, "TOCHE_H_hyper_mi5_corr_FINAL_LOST.rds")

set.seed(3)
TOCHE_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_TOCHE_M_A, name = "TOCHE_M")

saveRDS(TOCHE_M_hyper_mi, "TOCHE_M_hyper_mi5_corr_FINAL_LOST.rds")

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
names <- c("PC1", "PC2","PC3")

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

saveRDS(SAGU_H_hyper_mi, "SAGU_H_hyper_mi5_corr_FINAL_LOST.rds")

set.seed(3)
SAGU_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_SAGU_M_A, name = "SAGU_M")

saveRDS(SAGU_M_hyper_mi, "SAGU_M_hyper_mi5_corr_FINAL_LOST.rds")

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
names <- c("PC1", "PC2","PC3")

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

saveRDS(HON_H_hyper_mi, "HON_H_hyper_mi5_corr_FINAL_LOST.rds")

set.seed(3)
HON_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_HON_M_A, name = "HON_M")

saveRDS(HON_M_hyper_mi, "HON_M_hyper_mi5_corr_FINAL_LOST.rds")

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
names <- c("PC1", "PC2","PC3")

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

saveRDS(FLOR_H_hyper_mi, "FLOR_H_hyper_mi5_corr_FINAL_LOST.rds")

set.seed(3)
FLOR_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_FLOR_M_A, name = "FLOR_M")

saveRDS(FLOR_M_hyper_mi, "FLOR_M_hyper_mi5_corr_FINAL_LOST.rds")

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
names <- c("PC1", "PC2","PC3")

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

saveRDS(FUSA_H_hyper_mi, "FUSA_H_hyper_mi5_corr_FINAL_LOST.rds")

set.seed(3)
FUSA_M_hyper_mi <- hypervolume::hypervolume_svm(PCA_FUSA_M_A, name = "FUSA_M")

saveRDS(FUSA_M_hyper_mi, "FUSA_M_hyper_mi5_corr_FINAL_LOST.rds")

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
names <- c("PC1", "PC2","PC3")

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

TR <- read.csv("C:\\Users\\camil\\OneDrive - SELVA Investigación para la conservación en el neotropico\\Documentos\\Camila\\Publicaciones\\Chapman_Historical Occupancy\\All_species_traits_CorrJuli_FINAL_spLOST.csv", sep = ";")

TR <- TR %>% 
  filter(!HISTORICAL %in% c(0, NA))

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
dat <- read.csv("C:\\Users\\camil\\OneDrive - SELVA Investigación para la conservación en el neotropico\\Documentos\\Camila\\Publicaciones\\Chapman_Historical Occupancy\\SiteSpeciesMatrix_corr_LOST.csv", sep = ";",h = T)
head(dat)

#SAN ANTONIO

datANT <- subset(dat,Region == "SAN ANTONIO")

datANT <- datANT[order(dat$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species, Abundance, Period)

#Remove NAs 
datANT <- datANT[complete.cases(datANT), ]
datANT <- datANT[order(datANT$Species),]

tr_mi_z2 <- datANT[order(datANT$Species),]%>% 
  right_join(tr_mi_z1, by = "Species") %>% 
  select(Species,Abundance, Period,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z2 <- tr_mi_z2[complete.cases(tr_mi_z2), ]
tr_mi_z2 <- tr_mi_z2[order(tr_mi_z2$Species),] %>% 
  select(Species,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)

tr_mi_z2 <- tr_mi_z2[!duplicated(tr_mi_z2[,1]), ]

tr_mi_z2 <- data.frame(tr_mi_z2[,-1], row.names=tr_mi_z2[,1])


mat <- datANT %>% 
  pivot_wider(names_from =  Species, values_from = Abundance) %>% remove_rownames %>%
  column_to_rownames(var = "Period") 

mat[is.na(mat)] <- 0

#Estimating indices of functional diversity

fd_iNDICES <- dbFD(tr_mi_z2, mat,w.abun = FALSE, calc.FRic = TRUE, m = 7, stand.FRic = TRUE , calc.FDiv = TRUE,calc.CWM = FALSE, print.pco = FALSE, messages = TRUE)

fd_iNDICES



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
  select(Species,MODERN, HISTORICAL, ANT_H, ANT_M, ANT_EXT, POP_ANT,MIR_H_corr,MIR_M_LOST,MIR_EXT,POP_MIR,BARB_H_corr,BARB_M_LOST,BARB_EXT,POP_BARB,SAGU_H_corr, SAGU_M_LOST,SAGU_EXT,POP_SAGU, TOCHE_H_corr,TOCHE_M_LOST,TOCHE_EXT, POP_TOCHE, HON_H_corr,HON_M_LOST,HON_EXT, POP_HON, FLOR_H_corr, FLOR_M_LOST,FLOR_EXT, POP_FLOR,FUSA_H_corr, FUSA_M_LOST,FUSA_EXT, POP_FUSA,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) 

ANT <- trait_z %>% 
  filter(ANT_H >= 0)

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
## PCA ## SAN ANTONIO
#### --------------------------------------------------------------

#table just with columns of morphological traits and species names as rownames
tr_mi_z_ANT <- ANT %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z)%>% 
  column_to_rownames(var = "Species")

#Run PCA
pcaTotal_ANT <- princomp(tr_mi_z_ANT, cor = TRUE, scores = TRUE)
summary(pcaTotal_ANT)

#Scores
scoresTotal_ANT <- as.data.frame(pcaTotal_ANT$scores) %>% 
  tibble::rownames_to_column("Species") %>% 
  left_join(ANT, by = "Species")

#Community TPDc using PCAs as traits 

#Create new database for community TPD to compare extirpated, new additions and share specieswith Total PCA scores
TotalPCA<- scoresTotal_ANT %>% 
  mutate(SD = 1) %>%
  select(Species,Comp.1, Comp.2, Comp.3,POP_ANT,MODERN, HISTORICAL, SD) #POP = groups of species, A = Extirpated, B = Novel additions, C = Shared species

#PC1
TRA <- matrix(c(TotalPCA$Comp.1), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_ANT
ABUN <- TotalPCA %>% 
  select(Species, MODERN,POP_ANT) %>% 
  pivot_wider(names_from = Species, values_from = MODERN) %>% 
  column_to_rownames('POP_ANT')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC2
TRA <- matrix(c(TotalPCA$Comp.2), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_ANT
ABUN <- TotalPCA %>% 
  select(Species, HISTORICAL,POP_ANT) %>% 
  pivot_wider(names_from = Species, values_from = HISTORICAL) %>% 
  column_to_rownames('POP_ANT')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN )

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))


#PC3

TRA <- matrix(c(TotalPCA$Comp.3), ncol = 1)
SD <- matrix(c(TotalPCA$SD), ncol=1)
POP <- TotalPCA$POP_ANT
ABUN <- TotalPCA %>% 
  select(Species, HISTORICAL,POP_ANT) %>% 
  pivot_wider(names_from = Species, values_from = HISTORICAL) %>% 
  column_to_rownames('POP_ANT')

ABUN[is.na(ABUN)] <- 0

tpdmean <- TPDsMean(species = TotalPCA$Species, means = TRA, sds=SD, alpha = 1, samples = POP)


TPDc <- TPDc(TPDs = tpdmean, sampUnit = ABUN)

sapply(TPDc$TPDc$TPDc, sum)

plotTPD(TPD = TPDc, nRowCol = c(3,3))

#Density comparisons with kolmogorov-sANTnov test
library("sm")

A <- subset(TotalPCA, TotalPCA$POP == "A")
C <- subset(TotalPCA, TotalPCA$POP == "C")

#KS test
ks.test(A$Comp.1, C$Comp.1)# D = 0.34 P = 0.003

ks.test(A$Comp.2, C$Comp.2)# D = 0.18 P = 0.26

ks.test(A$Comp.3, C$Comp.3)# D = 0.24 P = 0.006


##Graphs of TPDs for PC1, 2 and 3

TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2) %>% 
  mutate(Comp3M = Comp.3)

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_ANT))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_ANT", values = c("red", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (65%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.8, 0.9))+
  scale_fill_manual(values=c("red",  "azure3"), 
                    name=NULL,breaks=c("A", "C"),labels=c("Extirpated species", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_ANT))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_ANT", values=c("red", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (18%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red", "azure3"), 
                    name=NULL,breaks=c("A", "C"),labels=c("Extirpated species", "Shared species"))

p

#PC3
windows()
p = ggplot(TotalPCA1, aes(x = Comp3M,fill = POP_ANT))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_ANT", values=c("red", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC3 - Bill length (6%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red", "azure3"), 
                    name=NULL,breaks=c("A", "C"),labels=c("Extirpated species", "Shared species"))

p



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

TotalPCA$POP_MIR[TotalPCA$POP_MIR == 'B'] <- 'C'

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
C <- subset(TotalPCA, TotalPCA$POP == "C")

#KS test
ks.test(A$Comp.1, C$Comp.1)# D = 32 P = 0.02

ks.test(A$Comp.2, C$Comp.2)# D = 0.16 P = 0.61

ks.test(A$Comp.3, C$Comp.3)# D = 0.26 P = 0.09


##Graphs of TPDs for PC1, 2 and 3

TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2) %>% 
  mutate(Comp3M = Comp.3)

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_MIR))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_MIR", values = c("red", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (59%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.8, 0.9))+
 scale_fill_manual(values=c("red",  "azure3"), 
                    name=NULL,breaks=c("A", "C"),labels=c("Extirpated species", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_MIR))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_MIR", values=c("red", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (22%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red", "azure3"), 
                    name=NULL,breaks=c("A", "C"),labels=c("Extirpated species", "Shared species"))

p

#PC3
windows()
p = ggplot(TotalPCA1, aes(x = Comp3M,fill = POP_MIR))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_MIR", values=c("red", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC3 - Bill length (6%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red", "azure3"), 
                    name=NULL,breaks=c("A", "C"),labels=c("Extirpated species", "Shared species"))

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

TotalPCA$POP_BARB[TotalPCA$POP_BARB == 'B'] <- 'C'

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

A <- subset(TotalPCA, TotalPCA$POP == "A")
C <- subset(TotalPCA, TotalPCA$POP == "C")

#KS test
ks.test(A$Comp.1, C$Comp.1)# D = 0.25 P = 0.06

ks.test(A$Comp.2, C$Comp.2)# D = 0.08 P = 0.97

ks.test(A$Comp.3, C$Comp.3)# D = 0.28 P = 0.02


##Graphs of TPDs for PC1, 2 and 3

TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2) %>% 
  mutate(Comp3M = Comp.3)

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_BARB))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_BARB", values = c("red", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (59%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.8, 0.9))+
  scale_fill_manual(values=c("red", "azure3"), 
                    name=NULL,breaks=c("A", "C"),labels=c("Extirpated species", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_BARB))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_BARB", values=c("red", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (22%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red", "azure3"), 
                    name=NULL,breaks=c("A", "C"),labels=c("Extirpated species", "Shared species"))

p

#PC3
windows()
p = ggplot(TotalPCA1, aes(x = Comp3M,fill = POP_BARB))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_BARB", values=c("red", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC3 - Bill length (7%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red", "azure3"), 
                    name=NULL,breaks=c("A", "C"),labels=c("Extirpated species", "Shared species"))

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

TotalPCA$POP_TOCHE[TotalPCA$POP_TOCHE == 'B'] <- 'C' 

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

A <- subset(TotalPCA, TotalPCA$POP == "A")
C <- subset(TotalPCA, TotalPCA$POP == "C")

#KS test
ks.test(A$Comp.1, C$Comp.1)#D = 0.28 P = 0.41

ks.test(A$Comp.2, C$Comp.2)#D = 0.32 P = 0.26

ks.test(A$Comp.3, C$Comp.3)#D = 0.23 P = 0.65


##Graphs of TPDs for PC1, 2 and 3

TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2) %>% 
  mutate(Comp3M = Comp.3)

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_TOCHE))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_TOCHE", values = c("red",  "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (65%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.75, 0.9))+ 
  scale_fill_manual(values=c("red",  "azure3"), 
                    name=NULL,breaks=c("A",  "C"),labels=c("Extirpated species",  "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_TOCHE))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_TOCHE", values=c("red",  "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (21%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red",  "azure3"), 
                    name=NULL,breaks=c("A", "C"),labels=c("Extirpated species", "Shared species"))

p

#PC3
windows()
p = ggplot(TotalPCA1, aes(x = Comp3M,fill = POP_TOCHE))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_TOCHE", values=c("red",  "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC3 - Bill length (6%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red",  "azure3"), 
                    name=NULL,breaks=c("A", "C"),labels=c("Extirpated species", "Shared species"))

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

TotalPCA$POP_HON[TotalPCA$POP_HON == 'B'] <- 'C' 

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

A <- subset(TotalPCA, TotalPCA$POP == "A")
C <- subset(TotalPCA, TotalPCA$POP == "C")

#KS test
ks.test(A$Comp.1, C$Comp.1)#D = 0.34 P = 0.17

ks.test(A$Comp.2, C$Comp.2)#D = 0.26 P = 0.43

ks.test(A$Comp.3, C$Comp.3)#D = 0.35 P = 0.14


##Graphs of TPDs for PC1, 2 and 3

TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2 ) %>% 
  mutate(Comp3M = Comp.3 )

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_HON))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_HON", values = c("red",  "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (63%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.75, 0.9))+ 
  scale_fill_manual(values=c("red",  "azure3"), 
                    name=NULL,breaks=c("A",  "C"),labels=c("Extirpated species", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_HON))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_HON", values=c("red", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (21%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red",  "azure3"), 
                    name=NULL,breaks=c("A",  "C"),labels=c("Extirpated species", "Shared species"))

p


#PC3
windows()
p = ggplot(TotalPCA1, aes(x = Comp3M,fill = POP_HON))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_HON", values=c("red", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC3 - Bill length (5%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red",  "azure3"), 
                    name=NULL,breaks=c("A",  "C"),labels=c("Extirpated species", "Shared species"))

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

TotalPCA$POP_FLOR[TotalPCA$POP_FLOR == 'B'] <- 'C'  

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

A <- subset(TotalPCA, TotalPCA$POP == "A")
C <- subset(TotalPCA, TotalPCA$POP == "C")

#KS test
ks.test(A$Comp.1, C$Comp.1)#D = 0.14 P = 0.20

ks.test(A$Comp.2, C$Comp.2)#D = 0.20 P = 0.03

ks.test(A$Comp.3, C$Comp.3)#D = 0.06 P = 0.98


##Graphs of TPDs for PC1, 2 and 3

TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2 ) %>% 
  mutate(Comp3M = Comp.3 )

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_FLOR))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_FLOR", values = c("red",  "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (64%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.75, 0.9))+ 
  scale_fill_manual(values=c("red",  "azure3"), 
                    name=NULL,breaks=c("A",  "C"),labels=c("Extirpated species", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_FLOR))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_FLOR", values=c("red",  "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (19%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red",  "azure3"), 
                    name=NULL,breaks=c("A",  "C"),labels=c("Extirpated species", "Shared species"))

p

#PC3
windows()
p = ggplot(TotalPCA1, aes(x = Comp3M,fill = POP_FLOR))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_FLOR", values=c("red",  "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC3 - Bill length (6%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red",  "azure3"), 
                    name=NULL,breaks=c("A",  "C"),labels=c("Extirpated species", "Shared species"))

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

TotalPCA$POP_SAGU[TotalPCA$POP_SAGU == 'B'] <- 'C'  

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

A <- subset(TotalPCA, TotalPCA$POP == "A")
C <- subset(TotalPCA, TotalPCA$POP == "C")

#KS test
ks.test(A$Comp.1, C$Comp.1)#D = 0.14 P = 0.21

ks.test(A$Comp.2, C$Comp.2)#D = 0.20 P = 0.02

ks.test(A$Comp.3, C$Comp.3)#D = 0.06 P = 0.98


##Graphs of TPDs for PC1, 2 and 3

TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2) %>% 
  mutate(Comp3M = Comp.3 )

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_SAGU))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_SAGU", values = c("red", "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (66%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.75, 0.9))+ 
  scale_fill_manual(values=c("red", "azure3"), 
                    name=NULL,breaks=c("A", "C"),labels=c("Extirpated species", "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_SAGU))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_SAGU", values=c("red",  "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (21%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red", "azure3"), 
                    name=NULL,breaks=c("A",  "C"),labels=c("Extirpated species",  "Shared species"))

p


#PC3
windows()
p = ggplot(TotalPCA1, aes(x = Comp3M,fill = POP_SAGU))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_SAGU", values=c("red",  "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC3 - Bill length (6%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+ 
  scale_fill_manual(values=c("red", "azure3"), 
                    name=NULL,breaks=c("A",  "C"),labels=c("Extirpated species",  "Shared species"))

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

TotalPCA$POP_FUSA[TotalPCA$POP_FUSA == 'B'] <- 'C'  

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

A <- subset(TotalPCA, TotalPCA$POP == "A")
C <- subset(TotalPCA, TotalPCA$POP == "C")

#KS test
ks.test(A$Comp.1, C$Comp.1)# D = 0.26P = 0.05

ks.test(A$Comp.2, C$Comp.2)# D = 0.17 P = 0.38

ks.test(A$Comp.3, C$Comp.3)# D = 0.13 P = 0.70


##Graphs of TPDs for PC1, 2 and 3

TotalPCA1 <- TotalPCA %>% 
  mutate(Comp2M = Comp.2 ) %>% 
  mutate(Comp3M = Comp.3 )

#PC1
windows()
p = ggplot(TotalPCA1, aes(x = Comp.1,fill = POP_FUSA))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_FUSA", values = c("red",  "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC1 - Body size (63%) ") + ylab("Trait prbability density") + ggtitle(" ")
p = p + theme(legend.position = c(0.75, 0.9))+ 
  scale_fill_manual(values=c("red",  "azure3"), 
                    name=NULL,breaks=c("A",  "C"),labels=c("Extirpated species",  "Shared species"))

p


#PC2
windows()
p = ggplot(TotalPCA1, aes(x = Comp2M,fill = POP_FUSA))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_FUSA", values=c("red",  "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC2 - Dispersal ability (22%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+
  scale_fill_manual(values=c("red",  "azure3"), 
                    name=NULL,breaks=c("A",  "C"),labels=c("Extirpated species",  "Shared species"))

p

#PC3
windows()
p = ggplot(TotalPCA1, aes(x = Comp3M,fill = POP_FUSA))
p = p + geom_density(alpha = 0.6) + scale_fill_manual(name = "POP_FUSA", values=c("red",  "azure3")) 
p = p + theme_bw(20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=20),axis.title=element_text(size=20), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))
p <- p + xlab("PC3 - Bill length (7%)") + ylab("Trait probability density") + ggtitle(" ")
p = p + theme(legend.position = "none")+
  scale_fill_manual(values=c("red",  "azure3"), 
                    name=NULL,breaks=c("A",  "C"),labels=c("Extirpated species",  "Shared species"))

p


##########################################
#Mapping vulnerable traits - extirpated species
################################################

#### Load data

# load: trait data
trait <- read.csv("C:\\Users\\camil\\OneDrive - SELVA Investigación para la conservación en el neotropico\\Documentos\\Camila\\Publicaciones\\Chapman_Historical Occupancy\\All_species_traits_CorrJuli_FINAL_spLOST.csv", sep = ";")

head(trait)
# Exploring trait correlations (verify trait column numbers)

trait <- trait %>% 
  filter(!HISTORICAL %in% c(0, NA))

# Exploring trait correlations (verify trait column numbers)

just_trait <- trait[, c(67,68,69,70,71,72,73,74,75,76,77)]
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
  select(Species,MODERN, HISTORICAL,ANT_EXT,MIR_EXT,BARB_EXT,SAGU_EXT,TOCHE_EXT, HON_EXT, FLOR_EXT, FUSA_EXT,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z, EXTINCT) 

ANT <- trait_z %>% 
  filter(ANT_EXT == 1 & 2)

MIR <- trait_z %>% 
  filter(MIR_EXT == 1 & 2)

BARB <- trait_z %>% 
  filter(BARB_EXT == 1 & 2)

SAGU <- trait_z %>% 
  filter(SAGU_EXT == 1 & 2)

TOCHE <- trait_z %>% 
  filter(TOCHE_EXT == 1 & 2)

HON <- trait_z %>% 
  filter(HON_EXT == 1 & 2)

FLOR <- trait_z %>% 
  filter(FLOR_EXT == 1 & 2)

FUSA <- trait_z %>% 
  filter(FUSA_EXT == 1 & 2)

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
dcc_TOTAL <- reshape2::melt(den_TOTAL$z)

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
  select(Species,ANT_EXT, MIR_EXT,BARB_EXT,SAGU_EXT,TOCHE_EXT, HON_EXT, FLOR_EXT, FUSA_EXT,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z, EXTINCT) %>% 
  left_join(scoresPCATotal_TOTAL, by = "Species") 

#### PCA plots
library(ggplot2)

# colour palette
col_pal <- colorRampPalette(c("grey30","grey60", "grey80", "white"))(200)


### Folowing lines indentify species that correspond to extinction,
#recolonization or new additions in different periods

extir_ANT <- subset(PCA_TOTAL, ANT_EXT == 1)
pextir_ANT <- subset(PCA_TOTAL, ANT_EXT == 2)

extir_MIR <- subset(PCA_TOTAL, MIR_EXT == 1)
pextir_MIR <- subset(PCA_TOTAL, MIR_EXT == 2)

extir_BARB <- subset(PCA_TOTAL, BARB_EXT == 1)
pextir_BARB <- subset(PCA_TOTAL, BARB_EXT == 2)

extir_SAGU <- subset(PCA_TOTAL, SAGU_EXT == 1)
pextir_SAGU <- subset(PCA_TOTAL, SAGU_EXT == 2)

extir_TOCHE <- subset(PCA_TOTAL, TOCHE_EXT == 1)
pextir_TOCHE <- subset(PCA_TOTAL, TOCHE_EXT == 2)

extir_HON <- subset(PCA_TOTAL, HON_EXT == 1)
pextir_HON <- subset(PCA_TOTAL, HON_EXT == 2)

extir_FLOR <- subset(PCA_TOTAL, FLOR_EXT == 1)
pextir_FLOR <- subset(PCA_TOTAL, FLOR_EXT == 2)

extir_FUSA <- subset(PCA_TOTAL, FUSA_EXT == 1)
pextir_FUSA <- subset(PCA_TOTAL, FUSA_EXT == 2)

# plot All Extirpated species
#PC1 and pc2
pca_plot_TOTAL <- ggplot(dcc_TOTAL, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_TOTAL, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_ANT, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_MIR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_BARB, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_SAGU, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_TOCHE, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_HON, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_FLOR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_FUSA, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  #geom_point(data = pextir_MIR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  #geom_point(data = pextir_BARB, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  #geom_point(data = pextir_SAGU, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  #geom_point(data = pextir_TOCHE, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  #geom_point(data = pextir_HON, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  #geom_point(data = pextir_FLOR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  #geom_point(data = pextir_FUSA, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
  
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
  labs(x = "PC1 - Size (64%)", y = "PC2 - Dispersal ability (20%)") +
  xlim(-8,18) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Extirpated species", color="black", size = 5)

# display plot
windows()
pca_plot_TOTAL


# plot Extirpated species - San antonio

pca_plot_TOTAL <- ggplot(dcc_TOTAL, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  geom_point(data = PCA_TOTAL, aes(x = Comp.1, y = Comp.2), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_ANT, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = pextir_ANT, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  labs(x = "PC1 - Size (64%)", y = "PC2 - Dispersal ability (20%)") +
  xlim(-8,15) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="San Antonio", color="black", size = 5)

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
  geom_point(data = pextir_MIR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  labs(x = "PC1 - Size (64%)", y = "PC2 - Dispersal ability (20%)") +
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
  geom_point(data = pextir_BARB, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  labs(x = "PC1 - Size (64%)", y = "PC2 - Dispersal ability (20%)") +
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
  geom_point(data = pextir_SAGU, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  labs(x = "PC1 - Size (64%)", y = "PC2 - Dispersal ability (20%)") +
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
  geom_point(data = pextir_TOCHE, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  labs(x = "PC1 - Size (64%)", y = "PC2 - Dispersal ability (20%)") +
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
  geom_point(data = pextir_HON, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  labs(x = "PC1 - Size (64%)", y = "PC2 - Dispersal ability (20%)") +
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
  geom_point(data = pextir_FLOR, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  labs(x = "PC1 - Size (64%)", y = "PC2 - Dispersal ability (20%)") +
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
  geom_point(data = pextir_FUSA, aes(x = Comp.1, y = Comp.2), size = 1.5, alpha = 0.8, colour = "orange") +
  
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
  labs(x = "PC1 - Size (64%)", y = "PC2 - Dispersal ability (20%)") +
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


#PC2 and pc3
#Subset for species in each Region and period

TOTAL_1 <- trait_z %>% 
  select(Species, weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z) %>% 
  right_join(scoresPCATotal_TOTAL, by = "Species") 

# kernel density estimation for each Region and each period
pc_raw_TOTAL <- TOTAL_1 %>% 
  # extract first two principal components
  dplyr::select(., Species, Comp.2, Comp.3) %>% 
  tibble::column_to_rownames(var = "Species")

# optimal bandwidth estimation
hpi_TOTAL <- Hpi(x = pc_raw_TOTAL)

# kernel density estimation Miraflores_Historical  
est_TOTAL <- kde(x = pc_raw_TOTAL, H = hpi_TOTAL, compute.cont = TRUE)  

den_TOTAL <- list(est_TOTAL$eval.points[[1]], est_TOTAL$eval.points[[2]], est_TOTAL$estimate)
names(den_TOTAL) <- c("x", "y", "z")
dimnames(den_TOTAL$z) <- list(den_TOTAL$x, den_TOTAL$y)
dcc_TOTAL <- reshape2::melt(den_TOTAL$z)

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
  select(Species,ANT_EXT,MIR_EXT,BARB_EXT,SAGU_EXT,TOCHE_EXT, HON_EXT, FLOR_EXT, FUSA_EXT,weight_z,T_cul_z, Bill_w_z, Wing_l_z,Tail_l_z, Tars_z, HWI_z, EXTINCT) %>% 
  left_join(scoresPCATotal_TOTAL, by = "Species") 

#### PCA plots
library(ggplot2)

# colour palette
col_pal <- colorRampPalette(c("grey30","grey60", "grey80", "white"))(200)


### Folowing lines indentify species that correspond to extinction,
#recolonization or new additions in different periods

extir_ANT <- subset(PCA_TOTAL, ANT_EXT == 1)
pextir_ANT <- subset(PCA_TOTAL, ANT_EXT == 2)

extir_MIR <- subset(PCA_TOTAL, MIR_EXT == 1)
pextir_MIR <- subset(PCA_TOTAL, MIR_EXT == 2)

extir_BARB <- subset(PCA_TOTAL, BARB_EXT == 1)
pextir_BARB <- subset(PCA_TOTAL, BARB_EXT == 2)

extir_SAGU <- subset(PCA_TOTAL, SAGU_EXT == 1)
pextir_SAGU <- subset(PCA_TOTAL, SAGU_EXT == 2)

extir_TOCHE <- subset(PCA_TOTAL, TOCHE_EXT == 1)
pextir_TOCHE <- subset(PCA_TOTAL, TOCHE_EXT == 2)

extir_HON <- subset(PCA_TOTAL, HON_EXT == 1)
pextir_HON <- subset(PCA_TOTAL, HON_EXT == 2)

extir_FLOR <- subset(PCA_TOTAL, FLOR_EXT == 1)
pextir_FLOR <- subset(PCA_TOTAL, FLOR_EXT == 2)

extir_FUSA <- subset(PCA_TOTAL, FUSA_EXT == 1)
pextir_FUSA <- subset(PCA_TOTAL, FUSA_EXT == 2)


pca_plot_TOTAL <- ggplot(dcc_TOTAL, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal)) +
  # points for species
  
  geom_point(data = PCA_TOTAL, aes(x = Comp.2, y = Comp.3), size = 0.3, alpha = 0.5, colour = "grey20") +
  geom_point(data = extir_MIR, aes(x = Comp.2, y = Comp.3), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_BARB, aes(x = Comp.2, y = Comp.3), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_SAGU, aes(x = Comp.2, y = Comp.3), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_TOCHE, aes(x = Comp.2, y = Comp.3), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_HON, aes(x = Comp.2, y = Comp.3), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_FLOR, aes(x = Comp.2, y = Comp.3), size = 1.5, alpha = 0.8, colour = "brown2") +
  geom_point(data = extir_FUSA, aes(x = Comp.2, y = Comp.3), size = 1.5, alpha = 0.8, colour = "brown2") +
  #geom_point(data = pextir_MIR, aes(x = Comp.2, y = Comp.3), size = 1.5, alpha = 0.8, colour = "orange") +
  #geom_point(data = pextir_BARB, aes(x = Comp.2, y = Comp.3), size = 1.5, alpha = 0.8, colour = "orange") +
  #geom_point(data = pextir_SAGU, aes(x = Comp.2, y = Comp.3), size = 1.5, alpha = 0.8, colour = "orange") +
  #geom_point(data = pextir_TOCHE, aes(x = Comp.2, y = Comp.3), size = 1.5, alpha = 0.8, colour = "orange") +
  #geom_point(data = pextir_HON, aes(x = Comp.2, y = Comp.3), size = 1.5, alpha = 0.8, colour = "orange") +
  #geom_point(data = pextir_FLOR, aes(x = Comp.2, y = Comp.3), size = 1.5, alpha = 0.8, colour = "orange") +
  #geom_point(data = pextir_FUSA, aes(x = Comp.2, y = Comp.3), size = 1.5, alpha = 0.8, colour = "orange") +
  
  
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.2, yend = Comp.3), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # probability kernels
  geom_contour(aes(z = value), breaks = cl_50_TOTAL, colour = "grey30", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_95_TOTAL, colour = "grey60", linewidth = 1) +
  geom_contour(aes(z = value), breaks = cl_99_TOTAL, colour = "grey70", linewidth = 1) +
  coord_equal() +
  # add arrows
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = Comp.2, yend = Comp.3), arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = loadings_sc_TOTAL, aes(x = 0, y = 0, xend = -Comp.2, yend = -Comp.3), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text_repel(data = loadings_sc_TOTAL, aes(x = Comp.2, y = Comp.3, label = trait), size = 4, nudge_x = 11, hjust = 0.5, direction = "y", segment.size = 0.5,segment.color = "grey89") +
  # axis labels - see comp_var
  labs(x = "PC2 - Dispersal ability (20%)", y = "PC3 - Bill morphology (6%)") +
  xlim(-8,18) +
  ylim(-8,8) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 22)) + 
  geom_text(x= 5, y=8, label="Extirpated species", color="black", size = 5)

# display plot
windows()
pca_plot_TOTAL

