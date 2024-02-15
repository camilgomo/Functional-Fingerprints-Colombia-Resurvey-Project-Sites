library(dplyr)
library(ggplot2)
library(AICcmodavg)
library(car)

dat <- read.csv("C:\\Users\\camil\\OneDrive - SELVA Investigación para la conservación en el neotropico\\Documentos\\Camila\\Publicaciones\\Chapman_Historical Occupancy\\glmm_table_Enero2024.csv", head = T, sep = ";")
head(dat)

#Exploring models

#Loss of natural vegetation 5

cand.models <- list()

cand.models[[1]] <- lm(dat$loss_percnat_5 ~ dat$extinct_spp)
cand.models[[2]]  <- lm(dat$loss_percnat_5 ~ dat$change_sp.richness)
cand.models[[3]]  <- lm(dat$loss_percnat_5 ~ dat$change_funrich)
cand.models[[4]]  <- lm(dat$loss_percnat_5 ~ dat$change_funeven)
cand.models[[5]]  <- lm(dat$loss_percnat_5 ~ dat$change_fundiver)
cand.models[[6]]  <- lm(dat$loss_percnat_5 ~ dat$change_fundisp)
cand.models[[7]]  <- lm(dat$loss_percnat_5 ~ dat$change_hypervol)
cand.models[[8]]  <- lm(dat$loss_percnat_5 ~ dat$change_hypervol_unique)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(cand.models), sep = " ")

##generate AICc table
aictab(cand.set = cand.models, modnames = Modnames, sort = TRUE)

##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = cand.models, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

summary(cand.models[[1]])
summary(cand.models[[2]])
summary(cand.models[[3]])
summary(cand.models[[4]])
summary(cand.models[[5]])
summary(cand.models[[6]])
summary(cand.models[[7]])
summary(cand.models[[8]])

#Loss of natural vegetation 20

cand.models2 <- list()

cand.models2[[1]] <- lm(dat$loss_percnat_20 ~ dat$extinct_spp)
cand.models2[[2]]  <- lm(dat$loss_percnat_20 ~ dat$change_sp.richness)
cand.models2[[3]]  <- lm(dat$loss_percnat_20 ~ dat$change_funrich)
cand.models2[[4]]  <- lm(dat$loss_percnat_20 ~ dat$change_funeven)
cand.models2[[5]]  <- lm(dat$loss_percnat_20 ~ dat$change_fundiver)
cand.models2[[6]]  <- lm(dat$loss_percnat_20 ~ dat$change_fundisp)
cand.models2[[7]]  <- lm(dat$loss_percnat_20 ~ dat$change_hypervol)
cand.models2[[8]]  <- lm(dat$loss_percnat_20 ~ dat$change_hypervol_unique)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(cand.models2), sep = " ")

##generate AICc table
aictab(cand.set = cand.models2, modnames = Modnames, sort = TRUE)

##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = cand.models2, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

summary(cand.models2[[1]])
summary(cand.models2[[2]])
summary(cand.models2[[3]])
summary(cand.models2[[4]])
summary(cand.models2[[5]])
summary(cand.models2[[6]])
summary(cand.models2[[7]])
summary(cand.models2[[8]])

#Loss of natural vegetation 50

cand.models3 <- list()

cand.models3[[1]] <- lm(dat$loss_percnat_50 ~ dat$extinct_spp)
cand.models3[[2]]  <- lm(dat$loss_percnat_50 ~ dat$change_sp.richness)
cand.models3[[3]]  <- lm(dat$loss_percnat_50 ~ dat$change_funrich)
cand.models3[[4]]  <- lm(dat$loss_percnat_50 ~ dat$change_funeven)
cand.models3[[5]]  <- lm(dat$loss_percnat_50 ~ dat$change_fundiver)
cand.models3[[6]]  <- lm(dat$loss_percnat_50 ~ dat$change_fundisp)
cand.models3[[7]]  <- lm(dat$loss_percnat_50 ~ dat$change_hypervol)
cand.models3[[8]]  <- lm(dat$loss_percnat_50 ~ dat$change_hypervol_unique)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(cand.models3), sep = " ")

##generate AICc table
aictab(cand.set = cand.models3, modnames = Modnames, sort = TRUE)

summary(cand.models3[[1]])
summary(cand.models3[[2]])
summary(cand.models3[[3]])
summary(cand.models3[[4]])
summary(cand.models3[[5]])
summary(cand.models3[[6]])
summary(cand.models3[[7]])
summary(cand.models3[[8]])

#Model Figures
library(ggrepel)
# New facet label names for scale variable
new.labs <- c("5 km", "20 km", "50 km")
names(new.labs) <- c("5", "20","50")

windows()
ggplot(dat, aes(x = extinct_spp, y = Loss_natVeg, group = factor(scale), color = factor(scale), label = Site )) + 
  geom_point(aes(color = factor(scale))) +
  geom_text_repel()+
  scale_color_manual(name = "Landscape scale", values = c("coral", "#5ab4ac","#d8b365"), labels =  c("5 km", "20 km","50 km"))+
  stat_smooth(method = "lm")+
  facet_wrap(~scale, labeller = labeller(scale = new.labs))+
  labs(x = "Extinct species", y = " % Loss of natural cover", title = "")+
  theme_bw(18) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

windows()
ggplot(dat, aes(x = change_sp.richness, y = Loss_natVeg, group = factor(scale), color = factor(scale), label = Site)) + 
  geom_point(aes(color = factor(scale))) +
  geom_text_repel()+
  scale_color_manual(name = "Landscape scale", values = c("coral", "#5ab4ac","#d8b365"), labels =  c("5 km", "20 km","50 km"))+
  stat_smooth(method = "lm")+
  facet_wrap(~scale, labeller = labeller(scale = new.labs))+
  labs(x = expression(paste(Delta," Species Richness")), y = "% Loss of natural cover", title = "")+
  theme_bw(18)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

windows()
ggplot(dat, aes(x = change_funrich, y = Loss_natVeg, group = factor(scale), color = factor(scale), label = Site)) + 
  geom_point(aes(color = factor(scale))) +
  geom_text_repel()+
  scale_color_manual(name = "Landscape scale", values = c("coral", "#5ab4ac","#d8b365"), labels =  c("5 km", "20 km","50 km"))+
  stat_smooth(method = "lm")+
  facet_wrap(~scale, labeller = labeller(scale = new.labs))+
  labs(x = expression(paste(Delta," Functional Richness")), y = "% Loss of natural cover", title = "")+
  theme_bw(18)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

windows()
ggplot(dat, aes(x = change_funeven, y = Loss_natVeg, group = factor(scale), color = factor(scale), label = Site)) + 
  geom_point(aes(color = factor(scale))) +
  geom_text_repel()+
  scale_color_manual(name = "Landscape scale", values = c("coral", "#5ab4ac","#d8b365"), labels =  c("5 km", "20 km","50 km"))+
  stat_smooth(method = "lm")+
  facet_wrap(~scale, labeller = labeller(scale = new.labs))+
  labs(x = expression(paste(Delta," Functional Evenness")), y = "% Loss of natural cover", title = "")+
  theme_bw(18)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

windows()
ggplot(dat, aes(x = change_fundiver, y = Loss_natVeg, group = factor(scale), color = factor(scale), label = Site)) + 
  geom_point(aes(color = factor(scale))) +
  geom_text_repel()+
  scale_color_manual(name = "Landscape scale", values = c("coral", "#5ab4ac","#d8b365"), labels =  c("5 km", "20 km","50 km"))+
  stat_smooth(method = "lm")+
  facet_wrap(~scale, labeller = labeller(scale = new.labs))+
  labs(x = expression(paste(Delta," Functional Divergence")), y = "% Loss of natural cover", title = "")+
  theme_bw(18)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

windows()
ggplot(dat, aes(x = change_fundisp, y = Loss_natVeg, group = factor(scale), color = factor(scale), label = Site)) + 
  geom_point(aes(color = factor(scale))) +
  geom_text_repel()+
  scale_color_manual(name = "Landscape scale", values = c("coral", "#5ab4ac","#d8b365"), labels =  c("5 km", "20 km","50 km"))+
  stat_smooth(method = "lm")+
  facet_wrap(~scale, labeller = labeller(scale = new.labs))+
  labs(x = expression(paste(Delta," Functional Dispersion")), y = "% Loss of natural cover", title = "")+
  theme_bw(18)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

windows()
ggplot(dat, aes(x = change_hypervol_unique, y = Loss_natVeg, group = factor(scale), color = factor(scale), label = Site)) + 
  geom_point(aes(color = factor(scale))) +
  geom_text_repel()+
  scale_color_manual(name = "Landscape scale", values = c("coral", "#5ab4ac","#d8b365"), labels =  c("5 km", "20 km","50 km"))+
  stat_smooth(method = "lm")+
  facet_wrap(~scale, labeller = labeller(scale = new.labs))+
  labs(x = expression(paste(Delta," Functional Unique Hypervolume")), y = "% Loss of natural cover", title = "")+
  theme_bw(18)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


windows()
ggplot(dat, aes(x = change_hypervol, y = Loss_natVeg, group = factor(scale), color = factor(scale), label = Site)) + 
  geom_point(aes(color = factor(scale))) +
  geom_text_repel()+
  scale_color_manual(name = "Landscape scale", values = c("coral", "#5ab4ac","#d8b365"), labels =  c("5 km", "20 km","50 km"))+
  stat_smooth(method = "lm")+
  facet_wrap(~scale, labeller = labeller(scale = new.labs))+
  labs(x = expression(paste(Delta," Functional Hypervolume")), y = "% Loss of natural cover", title = "")+
  theme_bw(18)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
