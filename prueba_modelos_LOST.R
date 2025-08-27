library(dplyr)
library(ggplot2)
library(AICcmodavg)
library(car)

dat <- read.csv("C:\\Users\\camil\\OneDrive - SELVA Investigación para la conservación en el neotropico\\Documentos\\Camila\\Publicaciones\\Chapman_Historical Occupancy\\glmm_table_Ago2025.csv", head = T, sep = ";")
head(dat)

#Summary visualization of changes
new.labs <- c("5 km", "20 km", "50 km")
names(new.labs) <- c("5", "20","50")
windows()
ggplot(dat, aes(x = Site, y = (change_funrich)*100)) + 
  geom_point(aes(size = 3)) +
  labs(x = "", y = expression(paste(Delta," Functional Richness", title = "")))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_text(x = 2, y = -100, label = "Historically higher")+
  geom_text(x = 2, y = 100, label = "Currently higher")+
  ylim(-100,100)+
  theme_bw(18) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5))

ggplot(dat, aes(x = Site, y = (change_funeven)*100)) + 
  geom_point(aes(size = 3)) +
  labs(x = "", y = expression(paste(Delta," Functional Evenness", title = "")))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_text(x = 2, y = -100, label = "Historically higher")+
  geom_text(x = 2, y = 100, label = "Currently higher")+
  ylim(-100,100)+
  theme_bw(18) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5))

ggplot(dat, aes(x = Site, y = (change_fundiver)*100)) + 
  geom_point(aes(size = 3)) +
  labs(x = "", y = expression(paste(Delta," Functional Divergence", title = "")))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_text(x = 2, y = -100, label = "Historically higher")+
  geom_text(x = 2, y = 100, label = "Currently higher")+
  ylim(-100,100)+
  theme_bw(18) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5))

ggplot(dat, aes(x = Site, y = (change_fundisp)*100)) + 
  geom_point(aes(size = 3)) +
  labs(x = "", y = expression(paste(Delta," Functional Dispersion", title = "")))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_text(x = 2, y = -100, label = "Historically higher")+
  geom_text(x = 2, y = 100, label = "Currently higher")+
  ylim(-100,100)+
  theme_bw(18) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5))

#Exploring models

#Percentage lost 1 (Extinct species)

cand.models <- list()

cand.models[[1]] <- lm(dat$percentage_lost_1 ~ dat$loss_percnat_5)
cand.models[[2]] <- lm(dat$percentage_lost_1 ~ dat$loss_percnat_20)
cand.models[[3]] <- lm(dat$percentage_lost_1 ~ dat$loss_percnat_50)

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

#Percentage lost 2 (Extinct species + possibly extinct spp)
cand.models2 <- list()

cand.models2[[1]] <- lm(dat$percentage_lost_2 ~ dat$loss_percnat_5)
cand.models2[[2]] <- lm(dat$percentage_lost_2 ~ dat$loss_percnat_20)
cand.models2[[3]] <- lm(dat$percentage_lost_2 ~ dat$loss_percnat_50)

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

# Change in Functional Richness
cand.models3 <- list()

cand.models3[[1]] <- lm((dat$change_funrich)*100 ~ dat$loss_percnat_5)
cand.models3[[2]] <- lm((dat$change_funrich)*100 ~ dat$loss_percnat_20)
cand.models3[[3]] <- lm((dat$change_funrich)*100 ~ dat$loss_percnat_50)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(cand.models3), sep = " ")

##generate AICc table
aictab(cand.set = cand.models3, modnames = Modnames, sort = TRUE)

##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = cand.models3, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

summary(cand.models3[[1]])
summary(cand.models3[[2]])
summary(cand.models3[[3]])

# Change in Functional Evenness
cand.models4 <- list()

cand.models4[[1]] <- lm((dat$change_funeven)*100 ~ dat$loss_percnat_5)
cand.models4[[2]] <- lm((dat$change_funeven)*100 ~ dat$loss_percnat_20)
cand.models4[[3]] <- lm((dat$change_funeven)*100 ~ dat$loss_percnat_50)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(cand.models4), sep = " ")

##generate AICc table
aictab(cand.set = cand.models4, modnames = Modnames, sort = TRUE)

##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = cand.models4, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

summary(cand.models4[[1]])
summary(cand.models4[[2]])
summary(cand.models4[[3]])

# Change in Functional Divergence
cand.models5 <- list()

cand.models5[[1]] <- lm((dat$change_fundiver)*100 ~ dat$loss_percnat_5)
cand.models5[[2]] <- lm((dat$change_fundiver)*100 ~ dat$loss_percnat_20)
cand.models5[[3]] <- lm((dat$change_fundiver)*100 ~ dat$loss_percnat_50)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(cand.models5), sep = " ")

##generate AICc table
aictab(cand.set = cand.models5, modnames = Modnames, sort = TRUE)

##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = cand.models5, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

summary(cand.models5[[1]])
summary(cand.models5[[2]])
summary(cand.models5[[3]])

# Change in Functional Dispersion
cand.models6 <- list()

cand.models6[[1]] <- lm((dat$change_fundisp)*100 ~ dat$loss_percnat_5)
cand.models6[[2]] <- lm((dat$change_fundisp)*100 ~ dat$loss_percnat_20)
cand.models6[[3]] <- lm((dat$change_fundisp)*100 ~ dat$loss_percnat_50)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(cand.models6), sep = " ")

##generate AICc table
aictab(cand.set = cand.models6, modnames = Modnames, sort = TRUE)

##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = cand.models6, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

summary(cand.models6[[1]])
summary(cand.models6[[2]])
summary(cand.models6[[3]])

# Change in Functional Hypervolume
cand.models7 <- list()

cand.models7[[1]] <- lm(dat$change_hypervol ~ dat$loss_percnat_5)
cand.models7[[2]] <- lm(dat$change_hypervol ~ dat$loss_percnat_20)
cand.models7[[3]] <- lm(dat$change_hypervol ~ dat$loss_percnat_50)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(cand.models7), sep = " ")

##generate AICc table
aictab(cand.set = cand.models7, modnames = Modnames, sort = TRUE)

##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = cand.models7, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

summary(cand.models7[[1]])
summary(cand.models7[[2]])
summary(cand.models7[[3]])


#Model Figures
library(ggrepel)
# New facet label names for scale variable
new.labs <- c("5 km", "20 km", "50 km")
names(new.labs) <- c("5", "20","50")

windows()
ggplot(dat, aes(x = Loss_natVeg, y = percentage_lost_1, group = factor(scale), color = factor(scale), label = Site )) + 
  geom_point(aes(color = factor(scale))) +
  geom_text_repel()+
  scale_color_manual(name = "Landscape scale", values = c("cyan4","darkorchid3", "red3"), labels =  c("5 km", "20 km","50 km"))+
  stat_smooth(method = "lm")+
  facet_wrap(~scale, labeller = labeller(scale = new.labs))+
  labs(y = "% Extirpated species", x = " % Loss of natural cover", title = "")+
  theme_bw(18) +
  ylim(-10,40)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

windows()
ggplot(dat, aes(y = percentage_lost_2, x = Loss_natVeg, group = factor(scale), color = factor(scale), label = Site )) + 
  geom_point(aes(color = factor(scale))) +
  geom_text_repel()+
  scale_color_manual(name = "Landscape scale", values = c("cyan4","darkorchid3", "red3"), labels =  c("5 km", "20 km","50 km"))+
  stat_smooth(method = "lm")+
  facet_wrap(~scale, labeller = labeller(scale = new.labs))+
  labs(y = "% Extirpated plus possibly extirpated species", x = " % Loss of natural cover", title = "")+
  theme_bw(18) +
  ylim(-10,40)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


windows()
ggplot(dat, aes(y = (change_funrich)*100, x = Loss_natVeg, group = factor(scale), color = factor(scale), label = Site)) + 
  geom_point(aes(color = factor(scale))) +
  geom_text_repel()+
  scale_color_manual(name = "Landscape scale", values = c("cyan4","darkorchid3", "red3"), labels =  c("5 km", "20 km","50 km"))+
  stat_smooth(method = "lm")+
  facet_wrap(~scale, labeller = labeller(scale = new.labs))+
  labs(y = expression(paste(Delta," Functional Richness")), x = "% Loss of natural cover", title = "")+
  theme_bw(18)+
  ylim(-100,100)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

windows()
ggplot(dat, aes(y = (change_funeven)*100, x = Loss_natVeg, group = factor(scale), color = factor(scale), label = Site)) + 
  geom_point(aes(color = factor(scale))) +
  geom_text_repel()+
  scale_color_manual(name = "Landscape scale", values = c("cyan4","darkorchid3", "red3"), labels =  c("5 km", "20 km","50 km"))+
  stat_smooth(method = "lm")+
  facet_wrap(~scale, labeller = labeller(scale = new.labs))+
  labs(y = expression(paste(Delta," Functional Evenness")), x = "% Loss of natural cover", title = "")+
  theme_bw(18)+
  ylim(-3,3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5))

windows()
ggplot(dat, aes(y = (change_fundiver)*100, x = Loss_natVeg, group = factor(scale), color = factor(scale), label = Site)) + 
  geom_point(aes(color = factor(scale))) +
  geom_text_repel()+
  scale_color_manual(name = "Landscape scale", values = c("cyan4","darkorchid3", "red3"), labels =  c("5 km", "20 km","50 km"))+
  stat_smooth(method = "lm")+
  facet_wrap(~scale, labeller = labeller(scale = new.labs))+
  labs(y = expression(paste(Delta," Functional Divergence")), x = "% Loss of natural cover", title = "")+
  theme_bw(18)+
  ylim(-3,3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        ,
        axis.text.x = element_text(angle = 45, vjust = 0.5))

windows()
ggplot(dat, aes(y = (change_fundisp)*100, x = Loss_natVeg, group = factor(scale), color = factor(scale), label = Site)) + 
  geom_point(aes(color = factor(scale))) +
  geom_text_repel()+
  scale_color_manual(name = "Landscape scale", values = c("cyan4","darkorchid3", "red3"), labels =  c("5 km", "20 km","50 km"))+
  stat_smooth(method = "lm")+
  facet_wrap(~scale, labeller = labeller(scale = new.labs))+
  labs(y = expression(paste(Delta," Functional Dispersion")), x = "% Loss of natural cover", title = "")+
  theme_bw(18)+
  ylim(-30,30)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

windows()
ggplot(dat, aes(y = change_hypervol, x = Loss_natVeg, group = factor(scale), color = factor(scale), label = Site)) + 
  geom_point(aes(color = factor(scale))) +
  geom_text_repel()+
  scale_color_manual(name = "Landscape scale", values = c("cyan4","darkorchid3", "red3"), labels =  c("5 km", "20 km","50 km"))+
  stat_smooth(method = "lm")+
  facet_wrap(~scale, labeller = labeller(scale = new.labs))+
  labs(y = expression(paste(Delta," Functional Hypervolume")), x = "% Loss of natural cover", title = "")+
  theme_bw(18)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
