library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(janitor)

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
  ylim(-2,2)+
  theme_bw(18) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5))

ggplot(dat, aes(x = Site, y = (change_fundiver)*100)) + 
  geom_point(aes(size = 3)) +
  labs(x = "", y = expression(paste(Delta," Functional Divergence", title = "")))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_text(x = 2, y = -100, label = "Historically higher")+
  geom_text(x = 2, y = 100, label = "Currently higher")+
  ylim(-2,2)+
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

cand1 <- lm(dat$percentage_lost_1 ~ dat$loss_percnat_5)
cand2 <- lm(dat$percentage_lost_1 ~ dat$loss_percnat_20)
cand3 <- lm(dat$percentage_lost_1 ~ dat$loss_percnat_50)

summary(cand1)
summary(cand2)
summary(cand3)

cand.predict1 <- cbind(dat,predict(cand1, interval = 'confidence'))
cand.predict2 <- cbind(dat,predict(cand2, interval = 'confidence'))
cand.predict3 <- cbind(dat,predict(cand3, interval = 'confidence'))


# plot the points (actual observations), regression line, and confidence interval
new.labs <- c("5 km", "20 km", "50 km")
names(new.labs) <- c("5", "20","50") 

p <- ggplot(cand.predict1, aes(loss_percnat_5, percentage_lost_1, label = Site))
p <- p + geom_point(aes(size = 2), color = "cyan4")
p <- p + geom_line(aes(loss_percnat_5, fit), color = "cyan4")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = "% Extirpated species", x = " % Loss of natural cover 5km", title = "")+
  theme_bw(18) +
  ylim(-10,40)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p


p <- ggplot(cand.predict2, aes(loss_percnat_20, percentage_lost_1, label = Site))
p <- p + geom_point(aes(size = 2), color = "darkorchid3")
p <- p + geom_line(aes(loss_percnat_20, fit), color = "darkorchid3")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = "% Extirpated species", x = " % Loss of natural cover 20km", title = "")+
  theme_bw(18) +
  ylim(-10,40)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p

p <- ggplot(cand.predict3, aes(loss_percnat_50, percentage_lost_1, label = Site))
p <- p + geom_point(aes(size = 2), color = "red3")
p <- p + geom_line(aes(loss_percnat_50, fit), color = "red3")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = "% Extirpated species", x = " % Loss of natural cover 50km", title = "")+
  theme_bw(18) +
  ylim(-10,40)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p

# Change in Functional Richness
cand4 <- lm((dat$change_funrich)*100 ~ dat$loss_percnat_5)
cand5 <- lm((dat$change_funrich)*100 ~ dat$loss_percnat_20)
cand6 <- lm((dat$change_funrich)*100 ~ dat$loss_percnat_50)

summary(cand4)
summary(cand5)
summary(cand6)

cand.predict4 <- cbind(dat,predict(cand4, interval = 'confidence'))
cand.predict5 <- cbind(dat,predict(cand5, interval = 'confidence'))
cand.predict6 <- cbind(dat,predict(cand6, interval = 'confidence'))


# plot the points (actual observations), regression line, and confidence interval

p <- ggplot(cand.predict4, aes(loss_percnat_5, (change_funrich)*100, label = Site))
p <- p + geom_point(aes(size = 2), color = "cyan4")
p <- p + geom_line(aes(loss_percnat_5, fit), color = "cyan4")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = expression(paste(Delta," Functional Richness")), x = "% Loss of natural cover 5km", title = "")+
  theme_bw(18) +
  ylim(-60,30)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p


p <- ggplot(cand.predict5, aes(loss_percnat_20, (change_funrich)*100, label = Site))
p <- p + geom_point(aes(size = 2), color = "darkorchid3")
p <- p + geom_line(aes(loss_percnat_20, fit), color = "darkorchid3")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = expression(paste(Delta," Functional Richness")), x = "% Loss of natural cover 20km", title = "")+
  theme_bw(18) +
  ylim(-60,30)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p

p <- ggplot(cand.predict6, aes(loss_percnat_50, (change_funrich)*100, label = Site))
p <- p + geom_point(aes(size = 2), color = "red3")
p <- p + geom_line(aes(loss_percnat_50, fit), color = "red3")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = expression(paste(Delta," Functional Richness")), x = "% Loss of natural cover 50km", title = "")+
  theme_bw(18) +
  ylim(-60,30)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
# Change in Functional Evenness

cand7 <- lm((dat$change_funeven)*100 ~ dat$loss_percnat_5)
cand8 <- lm((dat$change_funeven)*100 ~ dat$loss_percnat_20)
cand9 <- lm((dat$change_funeven)*100 ~ dat$loss_percnat_50)

summary(cand7)
summary(cand8)
summary(cand9)

cand.predict7 <- cbind(dat,predict(cand7, interval = 'confidence'))
cand.predict8 <- cbind(dat,predict(cand8, interval = 'confidence'))
cand.predict9 <- cbind(dat,predict(cand9, interval = 'confidence'))


# plot the points (actual observations), regression line, and confidence interval

p <- ggplot(cand.predict7, aes(loss_percnat_5, (change_funeven)*100, label = Site))
p <- p + geom_point(aes(size = 2), color = "cyan4")
p <- p + geom_line(aes(loss_percnat_5, fit), color = "cyan4")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = expression(paste(Delta," Functional Evenness")), x = "% Loss of natural cover 5km", title = "")+
  theme_bw(18) +
  ylim(-1.5,1.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p


p <- ggplot(cand.predict8, aes(loss_percnat_20, (change_funeven)*100, label = Site))
p <- p + geom_point(aes(size = 2), color = "darkorchid3")
p <- p + geom_line(aes(loss_percnat_20, fit), color = "darkorchid3")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = expression(paste(Delta," Functional Evenness")), x = "% Loss of natural cover 20km", title = "")+
  theme_bw(18) +
  ylim(-1.5,1.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p

p <- ggplot(cand.predict9, aes(loss_percnat_50, (change_funeven)*100, label = Site))
p <- p + geom_point(aes(size = 2), color = "red3")
p <- p + geom_line(aes(loss_percnat_50, fit), color = "red3")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = expression(paste(Delta," Functional Evenness")), x = "% Loss of natural cover 50km", title = "")+
  theme_bw(18) +
  ylim(-1.5,1.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p


# Change in Functional Divergence

cand10 <- lm((dat$change_fundiver)*100 ~ dat$loss_percnat_5)
cand11 <- lm((dat$change_fundiver)*100 ~ dat$loss_percnat_20)
cand12 <- lm((dat$change_fundiver)*100 ~ dat$loss_percnat_50)

summary(cand10)
summary(cand11)
summary(cand12)

cand.predict10 <- cbind(dat,predict(cand10, interval = 'confidence'))
cand.predict11 <- cbind(dat,predict(cand11, interval = 'confidence'))
cand.predict12 <- cbind(dat,predict(cand12, interval = 'confidence'))


# plot the points (actual observations), regression line, and confidence interval

p <- ggplot(cand.predict10, aes(loss_percnat_5, (change_fundiver)*100, label = Site))
p <- p + geom_point(aes(size = 2), color = "cyan4")
p <- p + geom_line(aes(loss_percnat_5, fit), color = "cyan4")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = expression(paste(Delta," Functional Divergence")), x = "% Loss of natural cover 5km", title = "")+
  theme_bw(18) +
  ylim(-2.5,2.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p


p <- ggplot(cand.predict11, aes(loss_percnat_20, (change_fundiver)*100, label = Site))
p <- p + geom_point(aes(size = 2), color = "darkorchid3")
p <- p + geom_line(aes(loss_percnat_20, fit), color = "darkorchid3")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = expression(paste(Delta," Functional Divergence")), x = "% Loss of natural cover 20km", title = "")+
  theme_bw(18) +
  ylim(-2.5,2.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p

p <- ggplot(cand.predict12, aes(loss_percnat_50, (change_fundiver)*100, label = Site))
p <- p + geom_point(aes(size = 2), color = "red3")
p <- p + geom_line(aes(loss_percnat_50, fit), color = "red3")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = expression(paste(Delta," Functional Divergence")), x = "% Loss of natural cover 50km", title = "")+
  theme_bw(18) +
  ylim(-2.5,2.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p


# Change in Functional Dispersion

cand13 <- lm((dat$change_fundisp)*100 ~ dat$loss_percnat_5)
cand14 <- lm((dat$change_fundisp)*100 ~ dat$loss_percnat_20)
cand15 <- lm((dat$change_fundisp)*100 ~ dat$loss_percnat_50)

summary(cand13)
summary(cand14)
summary(cand15)

cand.predict13 <- cbind(dat,predict(cand13, interval = 'confidence'))
cand.predict14 <- cbind(dat,predict(cand14, interval = 'confidence'))
cand.predict15 <- cbind(dat,predict(cand15, interval = 'confidence'))


# plot the points (actual observations), regression line, and confidence interval

p <- ggplot(cand.predict13, aes(loss_percnat_5, (change_fundisp)*100, label = Site))
p <- p + geom_point(aes(size = 2), color = "cyan4")
p <- p + geom_line(aes(loss_percnat_5, fit), color = "cyan4")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = expression(paste(Delta," Functional Dispersion")), x = "% Loss of natural cover 5km", title = "")+
  theme_bw(18) +
  ylim(-30,20)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p


p <- ggplot(cand.predict14, aes(loss_percnat_20, (change_fundisp)*100, label = Site))
p <- p + geom_point(aes(size = 2), color = "darkorchid3")
p <- p + geom_line(aes(loss_percnat_20, fit), color = "darkorchid3")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = expression(paste(Delta," Functional Dispersion")), x = "% Loss of natural cover 20km", title = "")+
  theme_bw(18) +
  ylim(-30,20)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p

p <- ggplot(cand.predict15, aes(loss_percnat_50, (change_fundisp)*100, label = Site))
p <- p + geom_point(aes(size = 2), color = "red3")
p <- p + geom_line(aes(loss_percnat_50, fit), color = "red3")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = expression(paste(Delta," Functional Dispersion")), x = "% Loss of natural cover 50km", title = "")+
  theme_bw(18) +
  ylim(-30,20)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p


# Change in Functional Hypervolume

cand16 <- lm(dat$change_hypervol ~ dat$loss_percnat_5)
cand17 <- lm(dat$change_hypervol ~ dat$loss_percnat_20)
cand18 <- lm(dat$change_hypervol ~ dat$loss_percnat_50)

summary(cand16)
summary(cand17)
summary(cand18)

cand.predict16 <- cbind(dat,predict(cand16, interval = 'confidence'))
cand.predict17 <- cbind(dat,predict(cand17, interval = 'confidence'))
cand.predict18 <- cbind(dat,predict(cand18, interval = 'confidence'))


# plot the points (actual observations), regression line, and confidence interval

p <- ggplot(cand.predict16, aes(loss_percnat_5, change_hypervol, label = Site))
p <- p + geom_point(aes(size = 2), color = "cyan4")
p <- p + geom_line(aes(loss_percnat_5, fit), color = "cyan4")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = expression(paste(Delta," Hypervolume")), x = "% Loss of natural cover 5km", title = "")+
  theme_bw(18) +
  ylim(-20,20)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p


p <- ggplot(cand.predict17, aes(loss_percnat_20, change_hypervol, label = Site))
p <- p + geom_point(aes(size = 2), color = "darkorchid3")
p <- p + geom_line(aes(loss_percnat_20, fit), color = "darkorchid3")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = expression(paste(Delta," Hypervolume")), x = "% Loss of natural cover 20km", title = "")+
  theme_bw(18) +
  ylim(-20,20)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p

p <- ggplot(cand.predict18, aes(loss_percnat_50, change_hypervol, label = Site))
p <- p + geom_point(aes(size = 2), color = "red3")
p <- p + geom_line(aes(loss_percnat_50, fit), color = "red3")
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) + geom_text_repel(stat = "unique")+
  labs(y = expression(paste(Delta," Hypervolume")), x = "% Loss of natural cover 50km", title = "")+
  theme_bw(18) +
  ylim(-30,20)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p

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
