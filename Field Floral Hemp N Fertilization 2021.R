######## R code by Mona Farnisa
######## Field Floral Hemp Data 2021
######## Floral hemp (Cannabis sativa L.) responses to nitrogen fertilization 
######## under field conditions in the high desert

#libraries to load
library(tidyverse)
library(broom)
library(lubridate)
library(dplyr)
library(car)
library(lme4) 
library(nlme)
library(glmmTMB)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(tidyverse)
library(multcompView)
library(emmeans)
library(multcomp)
library(see)
library(performance)
library(plyr)
library(multcompView)
library(emmeans)
library(multcomp)
library(gridExtra)
library(cowplot)
library(merTools)
library(broom)




############## Canopy Cover Field 2021 ###########
read.csv('Canopy Cover.csv')
field.canopy <- read.csv("Canopy Cover.csv", header = T, sep = ",")

#make dates into a date format  
field.canopy$date <- as.Date(field.canopy$Date,'%m/%d/%Y')

#create DAT column 
field.canopy$dat <- as.factor(difftime(as.POSIXct(field.canopy$date), as.POSIXct('2021-06-14', tz="PT"), units="days"))

#factor variety
field.canopy$variety <- as.factor(field.canopy$Variety)

#factor treatment 
field.canopy$treatment <- as.factor(field.canopy$Treatment)

#factor plots 
field.canopy$plot <- as.factor(field.canopy$Plot)

#factor rows 
field.canopy$row <- as.factor(field.canopy$Row)

# create blocks 
field.canopy <- field.canopy %>% 
  mutate(block = case_when(
    row == "1" ~ "1",
    row == "2" ~ "1", 
    row == "3" ~ "2", 
    row == "4" ~ "2", 
    row == "5" ~ "3",
    row == "6" ~ "3", 
    row == "7" ~ "4", 
    row == "8" ~ "4"
  ))

###cc per plant
field.canopy$ccperplant <- field.canopy$Canopy_m2.per.plant
field.canopy$vartreat <- as.factor(paste(field.canopy$variety, field.canopy$treatment))

###Check to see that all combinations are correct
xtabs(~ variety + date, data = field.canopy)
xtabs(~ treatment + dat, data = field.canopy)
xtabs(~ treatment + variety, data = field.canopy)
xtabs(~ row + variety, data = field.canopy)
xtabs(~ plot + variety, data = field.canopy)
xtabs(~ block + variety, data = field.canopy)

library(plyr)
#per plant
canopycover.stats <- ddply(field.canopy, c("vartreat", "dat", 'treatment', 'variety'),
                           summarise,
                           N    = length(ccperplant),#We use length instead of count.
                           mean = mean(ccperplant),
                           sd   = sd(ccperplant),
                           se   = sd / sqrt(N))
canopycover.stats
#write.table(canopycover.stats, file = 'Field CC means.csv', sep = ",", quote = FALSE, row.names = F)

levels(canopycover.stats$variety)[levels(canopycover.stats$variety)=='Tahoe'] <- 'Tahoe Cinco'

# presentation figure canopy cover 
ggplot(canopycover.stats, aes(x=dat, y=mean, group = vartreat, color = treatment, linetype = variety)) + 
  geom_point(size = 4) +
  geom_line(size = 1) +
  scale_linetype_manual(values=c('solid', "dotted", "dashed"))+
  scale_color_manual(values = c("#1B9E77", '#E7298A'))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, size=1,position=position_dodge(0)) +
  theme_classic() +
  ylab(bquote('Canopy cover  '(m^2~plant^-1))) +
  xlab('Days after transplanting') +
  labs(color = 'Treatment', linetype = 'Cultivar') + #rename legend title
  theme(axis.title.x = element_text(size = 16), axis.text = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  theme(legend.text=element_text(size=16)) + 
  #theme(legend.title = element_blank()) +
  theme(legend.title=element_text(size=14)) +
  theme(strip.text = element_text(size = 16)) +
  scale_y_continuous(expand = c(0,0),limits = c(0, 0.8))+
  #geom_text_repel(aes(label=l), size = 7, nudge_x = 0.2)+
  theme(legend.position = c(.1, .8)) + 
  theme(legend.box.background = element_rect(color="black", size=1))
  #theme(legend.position = 'none')

###### Model Building ~ canopy cover
c1 <- lmer(ccperplant ~ treatment*variety + (1|block) + (1|dat), data = field.canopy, REML = F)
summary(c1)
# check plots 
plot(c1)
par(mfrow=c(1,2))
qqnorm(resid(c1))
qqline(resid(c1))
# residuals
hist(residuals(c1), breaks=10)

c2 <- lmer(ccperplant ~ treatment*variety + (1|plot) + (1|dat), data = field.canopy, REML = F)
summary(c2)
# check plots 
plot(c2)
par(mfrow=c(1,2))
qqnorm(resid(c2))
qqline(resid(c2))
# residuals
hist(residuals(c2), breaks=10)

c3 <- lmer(ccperplant ~ treatment*variety + (1|block/plot) + (1|dat), data = field.canopy, REML = F)
summary(c3)
# check plots 
plot(c3)
par(mfrow=c(1,2))
qqnorm(resid(c3))
qqline(resid(c3))
# residuals
hist(residuals(c3), breaks=10)

c4 <- lmer(log(ccperplant) ~ treatment*variety + (1|block/plot) + (1|dat), data = field.canopy, REML = F)
summary(c4)
# check plots 
plot(c4)
par(mfrow=c(1,2))
qqnorm(resid(c4))
qqline(resid(c4))
# residuals
hist(residuals(c4), breaks=10)

c5 <- lmer(sqrt(ccperplant) ~ treatment*variety + (1|block/plot) + (1|dat), data = field.canopy, REML = F)
summary(c5)
# check plots 
plot(c5)
par(mfrow=c(1,2))
qqnorm(resid(c5))
qqline(resid(c5))
# residuals
hist(residuals(c5), breaks=10)

####### calculate slopes
coefficients(c5)

library(merTools)
RMSE.merMod(c1)
RMSE.merMod(c2)
RMSE.merMod(c3)
RMSE.merMod(c4)
RMSE.merMod(c5)

library(see)
library(performance)
### check assumptions
check_model(c5)

Anova(c5)

# estimated marginal means
e1 = emmeans(c5, specs = pairwise ~ treatment|variety, adjust = "none", type = "response") # compare levels of treatment within group
e1
e2 = emmeans(c5, ~ treatment*variety, adjust = "none")
e2
e3 = emmeans(c5, ~ treatment|variety, type = "response", adjust = "none")
e3
# comparisons 
# confidence intervals with statistical tests 
contrast(e2, method="pairwise",adjust="none", type = "response", infer = T)

pairs(e2, alpha = 0.05, adjust = "none") # contrast stats without CI's
#write.csv(SPAD.contrasts, file = 'SPAD.contrasts.csv', quote = FALSE, row.names = T)

#summary(contrast(e2, method="pairwise",adjust="tukey"), infer=c(TRUE, TRUE)) # same as pairs

pwpp(e2, by = "treatment", adjust = "none")
pairs(e2, adjust = "none")
plot(e2, comparisons = T, adjust = "none", by = "variety") # if red arrows overlap groups are not sig from each other

#cld(e2)
canopycover.lmer.cld <- cld(e2,
                            alpha   = 0.05,
                            reversed  = T, ### reverse the order of letters
                            Letters = letters,    ### Use lower-case letters for .group
                            #type    = "response", ### Report emmeans in orginal scale
                            adjust =  "none", ### no adjustment = LSD for multiple comparisons
                            method = "paiwise")    
canopycover.lmer.cld

# create labels 
canopycover.stats <- canopycover.stats %>% mutate(l = case_when(
  dat == "108" & vartreat == "Red Bordeaux Control" ~ "c", 
  dat == "108" & vartreat == "Berry Blossom N+" ~ "a",
  dat == "108" & vartreat == "Red Bordeaux N+" ~ "ab",
  dat == "108" & vartreat == "Tahoe Control" ~ "d",
  dat == "108" & vartreat == "Berry Blossom Control" ~ "c",
  dat == "108" & vartreat == "Tahoe N+" ~ "bc"))

######## black & white manuscript figure for canopy cover
cc.line.bw <- ggplot(canopycover.stats, aes(x=dat, y=mean, group = vartreat, shape = vartreat, linetype = vartreat)) + 
  geom_point(size = 4) +
  geom_line(size = 1) +
  scale_shape_manual(values=c(15, 0, 17, 2, 18, 9))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, size=1, position=position_dodge(0)) +
  theme_classic() +
  ylab(bquote('Canopy cover  '(m^2~plant^-1))) +
  xlab('Days after transplanting') +
  labs(shape = 'Cultivar - Treatment', linetype = 'Cultivar - Treatment') + #rename legend title
  theme(axis.title.x = element_text(size = 16), axis.text = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  #ggtitle('Canopy Cover - Field 2021') +
  theme(legend.text=element_text(size=16)) + 
  #theme(legend.title=element_text(size=14)) +
  theme(legend.title = element_blank()) +
  theme(strip.text = element_text(size = 16)) +
  scale_y_continuous(expand = c(0,0),limits = c(0, 0.8))+
  geom_text_repel(aes(label=l), size = 7, nudge_x = 0.2) +
  #theme(legend.position = c(.2, .8)) + 
  #theme(legend.box.background = element_rect(color="black", size=1)) +
  theme(legend.position = 'none')
cc.line.bw




############# Field Data 2021 Measurements (Plant Height, Stem Diameter) #####
read.csv('Field Data.csv')
measurements <- read.csv('Field Data.csv', header = T, sep = ",")

#make dates into a date format  
measurements$date <- as.Date(measurements$Date,'%m/%d/%Y')

#create DAT column 
measurements$dat <- as.factor(difftime(as.POSIXct(measurements$date), as.POSIXct('2021-06-14', tz="PT"), units="days"))

#factor variety
measurements$variety <- as.factor(measurements$Variety)

#factor treatment 
measurements$treatment <- as.factor(measurements$Treatment)

# factor row 
measurements$row <- as.factor(measurements$Row)

measurements$plant.number <- measurements$Plant.Number

# create blocks 
measurements <- measurements %>% 
  mutate(block = case_when(
    row == "1" ~ "1",
    row == "2" ~ "1", 
    row == "3" ~ "2", 
    row == "4" ~ "2", 
    row == "5" ~ "3",
    row == "6" ~ "3", 
    row == "7" ~ "4", 
    row == "8" ~ "4"
  ))

###Check to see that all combinations are correct
xtabs(~ variety + dat, data = measurements)
xtabs(~ treatment + dat, data = measurements)
xtabs(~ Plant.Number + dat, data = measurements)
xtabs(~ treatment + variety, data = measurements)
xtabs(~ block + variety, data = measurements)

measurements1 <- na.omit(measurements)

# remove days 44, 95
measurements1 <- subset(measurements1, dat != 44& dat != 72 & dat != 95 & dat != 109)

#Combined Variety & Treatment 
measurements1$vartreat <- as.factor(paste(measurements1$variety, measurements1$treatment))

height.stats <- ddply(measurements1, c("variety", 'treatment',"dat", 'vartreat'),
                      summarise,
                      N    = length(Height..cm.),#We use length instead of count.
                      mean = mean(Height..cm.),
                      sd   = sd(Height..cm.),
                      se   = sd / sqrt(N))
height.stats

library(ggplot2)
##Line graphs of mean and s.e.
ggplot(height.stats, aes(x=dat, y=mean, group = vartreat, color = treatment, linetype = variety)) + 
  geom_point(size = 4) +
  geom_line(size = 1) +
  scale_linetype_manual(values=c('solid', "dotted", "dashed"))+
  scale_color_manual(values = c("#1B9E77", '#E7298A'))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, size=1,position=position_dodge(0)) +
  theme_classic() +
  ylab(bquote('Plant height  '(cm~plant^-1))) +
  xlab('Days after transplanting') +
  labs(color = 'Treatment', linetype = 'Cultivar') + #rename legend title
  theme(axis.title.x = element_text(size = 16), axis.text = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  theme(legend.text=element_text(size=16)) + 
  #theme(legend.title = element_blank()) +
  theme(legend.title=element_text(size=14)) +
  theme(strip.text = element_text(size = 16)) +
  scale_y_continuous(expand = c(0,0),limits = c(0, 150))+
  #geom_text_repel(aes(label=l), size = 7, nudge_x = 0.2)+
  theme(legend.position = c(.1, .8)) + 
  theme(legend.box.background = element_rect(color="black", size=1))+ 
  theme(legend.position = 'none')

###### Model Building ~ Plant Height
h1 <- lmer(Height..cm. ~ treatment*variety + (1|block) + (1|dat), data = measurements1, REML = F)
h3 <- lmer(Height..cm. ~ treatment*variety + (1|block/Plant.Number) + (1|dat), data = measurements1, REML = F)
# model selection 
AIC(h1, h2, h3) 

h4 <- lmer(log(Height..cm.) ~ treatment*variety + (1|block/Plant.Number) + (1|dat), data = measurements1, REML = F)
h5 <- lmer(sqrt(Height..cm.) ~ treatment*variety + (1|block/Plant.Number) + (1|dat), data = measurements1, REML = F)

h4.1 <- lmer(log(Height..cm.) ~ treatment*variety + (1|block) + (1|dat), data = measurements1, REML = F)
h4.2 <- lmer(log(Height..cm.) ~ treatment*variety + (1|Plant.Number) + (1|dat), data = measurements1, REML = F)
h4.3 <- lmer(log(Height..cm.) ~ treatment+variety + (1|block/Plant.Number) + (1|dat), data = measurements1, REML = F)

library(merTools)
RMSE.merMod(h3)
RMSE.merMod(h4)
RMSE.merMod(h5)

library(see)
library(performance)
### check assumptions
#check_model(h3)
check_model(h4)
#check_model(h5)

Anova(h4)

# estimated marginal means
e1 = emmeans(h4, specs = pairwise ~ treatment|variety, adjust = "none", type = "response") # compare levels of treatment within group
e1
e2 = emmeans(h4, ~ treatment*variety, adjust = "none")
e2
e3 = emmeans(h4, ~ treatment|variety, type = "response", adjust = "none")
e3
# comparisons 
# confidence intervals with statistical tests 
contrast(e2, method="pairwise",adjust="none", type = "response", infer = T)

SPAD.contrasts <- pairs(e2, alpha = 0.05, adjust = "none") # contrast stats without CI's
#write.csv(SPAD.contrasts, file = 'SPAD.contrasts.csv', quote = FALSE, row.names = T)
#summary(contrast(e2, method="pairwise",adjust="tukey"), infer=c(TRUE, TRUE)) # same as pairs

pwpp(e2, by = "treatment", adjust = "none")
pairs(e2, adjust = "none")
plot(e2, comparisons = T, adjust = "none", by = "variety") # if red arrows overlap groups are not sig from each other

cld(e2,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")    ### Tukey adjustment for multiple comparisons

# create labels 
height.stats <- height.stats %>% mutate(l = case_when(
  dat == "102" & vartreat == "Red Bordeaux Control" ~ "c", 
  dat == "102" & vartreat == "Berry Blossom N+" ~ "a",
  dat == "102" & vartreat == "Red Bordeaux N+" ~ "b",
  dat == "102" & vartreat == "Tahoe  Control" ~ "d",
  dat == "102" & vartreat == "Berry Blossom Control" ~ "cd",
  dat == "102" & vartreat == "Tahoe  N+" ~ "ab"))

########### manuscript plant height figure 
line.height.bw <- ggplot(height.stats, aes(x=dat, y=mean, group = vartreat, shape = vartreat, linetype = vartreat)) + 
  geom_point(size = 4) +
  geom_line(size = 1) +
  scale_shape_manual(labels = c('BBc', 'BBn', 'RBc', 'RBn', 'TCc', 'TCn'), values=c(15, 0, 17, 2, 18, 9))+
  scale_linetype_discrete(labels = c('BBc', 'BBn', 'RBc', 'RBn', 'TCc', 'TCn')) +
  #geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, size=1,position=position_dodge(0)) +
  theme_classic() +
  scale_fill_grey()+
  ylab(bquote('Plant height  '(cm~plant^-1))) +
  xlab('Days after transplanting') +
  #labs(shape = 'Cultivar - Treatment', linetype = 'Cultivar - Treatment') + #rename legend title
  theme(axis.title.x = element_text(size = 16), axis.text = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  #ggtitle('Plant Height linear model - Field 2021') +
  theme(legend.text=element_text(size=16)) + 
  theme(legend.title = element_blank()) +
  #theme(legend.title=element_text(size=14)) +
  theme(strip.text = element_text(size = 16)) +
  scale_y_continuous(expand = c(0,0),limits = c(0, 150))+
  geom_text_repel(aes(label=l), size = 7, nudge_x = 0.2)+
  theme(legend.position = c(.1, .8)) + 
  theme(legend.box.background = element_rect(color="black", size=1))
#theme(legend.position = 'none')
line.height.bw

###### merge canopy cover and plant height manuscript figures ####### 
cc.height.bw2 <- plot_grid(cc.line.bw, line.height.bw, labels = "AUTO", label_size = 18)

#ggsave(cc.height.bw2, file="cc.height.bw.TIFF",width=13, height=6,dpi=600,path="C:/Users/mfarnisa/OneDrive - University of Nevada, Reno/Thesis/Figures")


######### Model Building - Stem Diameter Stats
diameter.stats <- ddply(measurements1, c("variety", "treatment", "dat", "vartreat"),
                        summarise,
                        N    = length(Stem.Diameter..mm.),#We use length instead of count.
                        mean = mean(Stem.Diameter..mm.),
                        sd   = sd(Stem.Diameter..mm.),
                        se   = sd / sqrt(N))
diameter.stats
#write.table(diameter.stats, file = 'Stem Diameter Means field 2021.csv', sep = ",", quote = FALSE, row.names = F)

ggplot(diameter.stats, aes(x=dat, y=mean, group = vartreat, shape = vartreat, linetype = vartreat)) + 
  geom_point(size = 4) +
  geom_line(size = 1) +
  scale_shape_manual(labels = c('BBc', 'BBn', 'RBc', 'RBn', 'Tc', 'Tn'), values=c(15, 0, 17, 2, 18, 9))+
  scale_linetype_discrete(labels = c('BBc', 'BBn', 'RBc', 'RBn', 'Tc', 'Tn')) +
  #geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, size=1,position=position_dodge(0)) +
  theme_classic() +
  scale_fill_grey()+
  ylab(bquote('Diameter  '(mm~plant^-1))) +
  xlab('Days after transplanting') +
  #labs(shape = 'Cultivar - Treatment', linetype = 'Cultivar - Treatment') + #rename legend title
  theme(axis.title.x = element_text(size = 16), axis.text = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  #ggtitle('Plant Height linear model - Field 2021') +
  theme(legend.text=element_text(size=16)) + 
  theme(legend.title = element_blank()) +
  #theme(legend.title=element_text(size=14)) +
  theme(strip.text = element_text(size = 16)) +
  #scale_y_continuous(expand = c(0,0),limits = c(0, 150))+
  #geom_text_repel(aes(label=l), size = 7, nudge_x = 0.2)+
  theme(legend.position = c(.1, .8)) + 
  theme(legend.box.background = element_rect(color="black", size=1))
#theme(legend.position = 'none')

library(lme4)
library(rcompanion)
d.transfrom = transformTukey(measurements1$Stem.Diameter..mm., plotit=T)
d.sqr = (measurements1$Stem.Diameter..mm.)^(0.4)   # Avoid complex numbers

d1 <- lmer(Stem.Diameter..mm. ~ treatment*variety + (1|block) + (1|dat), data = measurements1, REML = F)
#d3 <- lmer(Stem.Diameter..mm. ~ treatment*variety + (1|block/Plant.Number) + (1|dat), data = measurements1, REML = F)
d4 <- lmer(log(Stem.Diameter..mm.) ~ treatment*variety + (1|block) + (1|dat), data = measurements1, REML = F)
d5 <- lmer(d.sqr ~ treatment*variety + (1|block) + (1|dat), data = measurements1, REML = F)

library(merTools)
RMSE.merMod(d1)
RMSE.merMod(d4)
RMSE.merMod(d5)

library(see)
library(performance)
### diagnostic plots
#check_model(d1)
check_model(d4)
check_model(d5)

#Anova(d1)
Anova(d4)
Anova(d5)

# estimated marginal means
e1 = emmeans(d4, specs = pairwise ~ treatment|variety, adjust = "none", type = "response") # compare levels of treatment within group
e1
e2 = emmeans(d4, ~ treatment*variety, type = "response", adjust = "none")
e2
e3 = emmeans(d5, ~ treatment*variety, type = "response", adjust = "none")
e3
# comparisons 
# confidence intervals with statistical tests 
contrast(e2, method="pairwise",adjust="none", type = "response", infer = T)

SPAD.contrasts <- pairs(e2, alpha = 0.05, adjust = "none") # contrast stats without CI's
#write.csv(SPAD.contrasts, file = 'SPAD.contrasts.csv', quote = FALSE, row.names = T)
#summary(contrast(e2, method="pairwise",adjust="tukey"), infer=c(TRUE, TRUE)) # same as pairs

pwpp(e2, by = "treatment", adjust = "none")
pairs(e2, adjust = "none")
plot(e2, comparisons = T, adjust = "none", by = "variety") # if red arrows overlap groups are not sig from each other

cld(e2,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    #type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")    ### no adjustment = LSD for multiple comparisons
cld(e3,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")    ### no adjustment = LSD for multiple comparisons



############# Dried Shoot Biomass ##############
read.csv('Biomass Weights.csv')
biomass <- read.csv('Biomass Weights.csv', header = T, sep = ",")

#make dates into a date format  
biomass$date <- as.Date(biomass$Date,'%m/%d/%Y')

#create DAT column 
biomass$dat <- as.factor(difftime(as.POSIXct(biomass$date), as.POSIXct('2021-06-14', tz="PT"), units="days"))

#factor variety
biomass$variety <- as.factor(biomass$Variety)

#factor treatment 
biomass$treatment <- as.factor(biomass$Treatment)

# factor row
biomass$row <- as.factor(biomass$Row)

# create blocks 
biomass <- biomass %>% 
  mutate(block = case_when(
    row == "1" ~ "1",
    row == "2" ~ "1", 
    row == "3" ~ "2", 
    row == "4" ~ "2", 
    row == "5" ~ "3",
    row == "6" ~ "3", 
    row == "7" ~ "4", 
    row == "8" ~ "4"
  ))

biomass$vartreat <- as.factor(paste(biomass$variety, biomass$treatment))

###Check to see that all combinations are correct
xtabs(~ variety + date, data = biomass)
xtabs(~ treatment + dat, data = biomass)
xtabs(~ treatment + variety, data = biomass)
xtabs(~ block + variety, data = biomass)

#calculate summary statistics of variables 
leafshoot.stats <- ddply(biomass, c("variety", "treatment"),
                         summarise,
                         N    = length(Leaf.Shoot.Total.Dried_g),#We use length instead of count.
                         mean = mean(Leaf.Shoot.Total.Dried_g),
                         sd   = sd(Leaf.Shoot.Total.Dried_g),
                         se   = sd / sqrt(N))
leafshoot.stats
#write.table(leafshoot.stats, file = 'Field Shoot Biomass Means.csv', sep = ",", quote = FALSE, row.names = F)
levels(leafshoot.stats$variety)[levels(leafshoot.stats$variety)=='Tahoe '] <- 'Tahoe Cinco'

boxplot(Leaf.Shoot.Total.Dried_g ~ vartreat, data = biomass, col = "lightgray")

library(ggplot2)
#Bar plot of mean and s.e.
biomass.treatment.bar <- ggplot(leafshoot.stats, aes(x = treatment, y=mean))+
  facet_wrap(~ variety) +
  #facet_grid(vars(Leaf.Position), vars(variety), scales = "free", space = 'free') + 
  geom_bar(stat = 'identity', position="dodge")+ 
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab("Shoot Weight (g)") + 
  xlab("Treatment") +
  ggtitle('Dried Biomass - Field 2021') + 
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title=element_text(size=14)) + 
  theme(strip.text = element_text(size = 14)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1500)) + 
  #scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = expansion(add=c(0.2,0.2)))
#geom_text(aes(label=letters.biomass, y = mean + se), vjust = -.5, size = 6)
biomass.treatment.bar

###### Model Building ~ biomass #######
b1 <- lmer(Leaf.Shoot.Total.Dried_g ~ treatment*variety + (1|block), data = biomass, REML = F)
# b2 <- lmer(log(Leaf.Shoot.Total.Dried_g) ~ treatment*variety + (1|block), data = biomass, REML = F)
# b3 <- lmer(sqrt(Leaf.Shoot.Total.Dried_g) ~ treatment*variety + (1|block), data = biomass, REML = F)

b4 <- lmer(Leaf.Shoot.Total.Dried_g ~ treatment+variety + (1|block), data = biomass, REML = F)
# b5 <- lmer(log(Leaf.Shoot.Total.Dried_g) ~ treatment+variety + (1|block), data = biomass, REML = F)
# b6 <- lmer(sqrt(Leaf.Shoot.Total.Dried_g) ~ treatment+variety + (1|block), data = biomass, REML = F)

b7 <- lmer(leaf.sqr ~ treatment*variety + (1|block), data = biomass, REML = F)

interaction.plot(x.factor     = biomass$treatment,
                 trace.factor = biomass$variety,
                 response     = biomass$Leaf.Shoot.Total.Dried_g,
                 fun = mean,
                 type="b",
                 col=c("black","red","green"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

library(rcompanion)
T_tuk = transformTukey(biomass$Leaf.Shoot.Total.Dried_g, plotit=T)

leaf.sqr = (biomass$Leaf.Shoot.Total.Dried_g)^(0.35)   # Avoid complex numbers

### check assumptions
check_model(b1)
#check_model(b4)
check_model(b7)

AIC(b1, b4)

Anova(b1)
#Anova(b4)
Anova(b7)

# estimated marginal means
e2 = emmeans(b1, ~ treatment*variety, type = "response", adjust = "none")
e2
e3 = emmeans(b7, ~ treatment*variety, type = "response", adjust = "none")
e3

cld(e2,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    #type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")    ### no adjustment = LSD for multiple comparisons
cld(e3,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    #type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")    ### no adjustment = LSD for multiple comparisons

# comparisons 
# confidence intervals with statistical tests 
contrast(e3, method="pairwise",adjust="none", type = "response", infer = T)

pairs(e2, alpha = 0.05, adjust = "none") # contrast stats without CI's
#write.csv(SPAD.contrasts, file = 'SPAD.contrasts.csv', quote = FALSE, row.names = T)

#summary(contrast(e2, method="pairwise",adjust="tukey"), infer=c(TRUE, TRUE)) # same as pairs

pwpp(e2, adjust = "none")
pairs(e2, adjust = "none")
plot(e2, comparisons = T, adjust = "none", by = "variety") # if red arrows overlap groups are not sig from each other


letters.biomass = c("c", "a", "c", "b", "d", "bc")

ggplot(leafshoot.stats, aes(x=treatment, y=mean, fill = treatment)) + 
  facet_wrap(~ variety) +
  geom_bar(stat = 'identity', position="dodge")+ 
  scale_fill_manual(values = c("#1B9E77", '#E7298A'))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, size=1,position=position_dodge(0)) +
  theme_classic() +
  ylab(bquote('Dried shoot biomass  '(g~plant^-1))) +
  xlab('Treatment') +
  labs(color = 'Treatment', linetype = 'Cultivar') + #rename legend title
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title=element_text(size=14)) + 
  theme(strip.text = element_text(size = 14)) + 
  theme(axis.title.y = element_text(vjust=0.9)) +
  scale_y_continuous(expand = c(0,0),limits = c(0, 1600))+
  scale_x_discrete(expand = expansion(add=c(0.2,0.2)))+
  geom_text(aes(label=letters.biomass, y = mean + se), vjust = -.5, size = 6)+
  theme(legend.position = 'none')



############# Dried Inflorescence Biomass ################

#calculate summary statistics of variables 
flower.stats <- ddply(biomass, c("variety", "treatment"),
                      summarise,
                      N    = length(Flower.Total.Dried_g),#We use length instead of count.
                      mean = mean(Flower.Total.Dried_g),
                      sd   = sd(Flower.Total.Dried_g),
                      se   = sd / sqrt(N))
flower.stats
#write.table(flower.stats, file = 'Field Flower Biomass Means.csv', sep = ",", quote = FALSE, row.names = F)

levels(flower.stats$variety)[levels(flower.stats$variety)=='Tahoe '] <- 'Tahoe Cinco'

flower.treatment.bar <- ggplot(flower.stats, aes(x = treatment, y=mean))+
  facet_wrap(~ variety) +
  #facet_grid(vars(Leaf.Position), vars(variety), scales = "free", space = 'free') + 
  geom_bar(stat = 'identity', position="dodge")+ 
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab("Flower Weight (g)") + 
  xlab("Treatment") +
  ggtitle('Dried Biomass - Field 2021') + 
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title=element_text(size=14)) + 
  theme(strip.text = element_text(size = 14)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,900)) + 
  #scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = expansion(add=c(0.2,0.2)))
  #geom_text(aes(label=letters.flower, y = mean + se), vjust = -.5, size = 6)
flower.treatment.bar

interaction.plot(x.factor     = biomass$treatment,
                 trace.factor = biomass$variety,
                 response     = biomass$Flower.Total.Dried_g,
                 fun = mean,
                 type="b",
                 col=c("black","red","green"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

f1 <- lmer(Flower.Total.Dried_g ~ treatment*variety + (1|block), data = biomass, REML = F)
f4 <- lmer(flower.sqr ~ treatment*variety + (1|block), data = biomass, REML = F)

#f2 <- lmer(log(Flower.Total.Dried_g) ~ treatment*variety + (1|block), data = biomass, REML = F)
#f3 <- lmer(sqrt(Flower.Total.Dried_g) ~ treatment*variety + (1|block), data = biomass, REML = F)

library(rcompanion)
T_tuk = transformTukey(biomass$Flower.Total.Dried_g, plotit=T)

flower.sqr = (biomass$Flower.Total.Dried_g)^(0.25)   # Avoid complex numbers

### check assumptions
check_model(f1)
#check_model(f2)
#check_model(f3)
check_model(f4)

Anova(f1)
#Anova(f2)
#Anova(f3)
Anova(f4)

# estimated marginal means
ff = emmeans(f1, ~ treatment*variety, type = "response", adjust = "none")
ff
ff3 = emmeans(f4, ~ treatment*variety, type = "response", adjust = "none")
ff3
# comparisons 
# confidence intervals with statistical tests 
contrast(ff, method="pairwise",adjust="none", type = "response", infer = T)

pairs(e2, alpha = 0.05, adjust = "none") # contrast stats without CI's
#write.csv(SPAD.contrasts, file = 'SPAD.contrasts.csv', quote = FALSE, row.names = T)

#summary(contrast(e2, method="pairwise",adjust="tukey"), infer=c(TRUE, TRUE)) # same as pairs

pwpp(e2, adjust = "none")
pairs(e2, adjust = "none")
plot(e2, comparisons = T, adjust = "none", by = "variety") # if red arrows overlap groups are not sig from each other

cld(ff,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    #type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")    ### no adjustment = LSD for multiple comparisons
cld(ff3,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    #type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")    ### no adjustment = LSD for multiple comparisons

#write.csv(SPAD.lmer.cld, file = 'SPAD.cld.csv', quote = FALSE, row.names = T)
letters.flower = c("cd", "ab", "cd", "bc", "d", "a")

ggplot(flower.stats, aes(x=treatment, y=mean, fill = treatment)) + 
  facet_wrap(~ variety) +
  geom_bar(stat = 'identity', position="dodge")+ 
  scale_fill_manual(values = c("#1B9E77", '#E7298A'))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, size=1,position=position_dodge(0)) +
  theme_classic() +
  ylab(bquote('Dried inflorescence biomass  '(g~plant^-1))) +
  xlab('Treatment') +
  labs(color = 'Treatment', linetype = 'Cultivar') + #rename legend title
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title=element_text(size=14)) + 
  theme(strip.text = element_text(size = 14)) + 
  theme(axis.title.y = element_text(vjust=0.9)) +
  scale_y_continuous(expand = c(0,0),limits = c(0, 1000))+
  scale_x_discrete(expand = expansion(add=c(0.2,0.2)))+
  geom_text(aes(label=letters.flower, y = mean + se), vjust = -.5, size = 6)+
  theme(legend.position = 'none')

########### Dried Inflorescence-to-Shoot Biomass Ratio ###########
library(MASS)
biomass$ratio = biomass$Flower.Total.Dried_g/biomass$Leaf.Shoot.Total.Dried_g

#calculate summary statistics of variables 
ratio.stats <- ddply(biomass, c("variety", "treatment"),
                     summarise,
                     N    = length(ratio),#We use length instead of count.
                     mean = mean(ratio),
                     sd   = sd(ratio),
                     se   = sd / sqrt(N))
ratio.stats
#write.table(ratio.stats, file = 'Biomass Ratio Means.csv', sep = ",", quote = FALSE, row.names = F)

ggplot(ratio.stats, aes(x=treatment, y=mean, fill = treatment)) + 
  facet_wrap(~ variety) +
  geom_point(stat = 'identity', position="dodge")+ 
  scale_fill_manual(values = c("#1B9E77", '#E7298A'))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, size=1,position=position_dodge(0)) +
  theme_classic() +
  ylab(bquote('Dried inflorescence biomass  '(g~plant^-1))) +
  xlab('Treatment') +
  labs(color = 'Treatment', linetype = 'Cultivar') + #rename legend title
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title=element_text(size=14)) + 
  theme(strip.text = element_text(size = 14)) + 
  theme(axis.title.y = element_text(vjust=0.9)) +
  scale_y_continuous(expand = c(0,0),limits = c(0, 2))+
  scale_x_discrete(expand = expansion(add=c(0.2,0.2)))+
  #geom_text(aes(label=letters.biomass, y = mean + se), vjust = -.5, size = 6)+
  theme(legend.position = 'none')

r1 <- lmer(ratio ~ treatment*variety + (1|row), data = biomass, REML = F)
r2 <- lmer(log(ratio) ~ treatment*variety + (1|row), data = biomass, REML = F)
r3 <- lmer(sqrt(ratio) ~ treatment*variety + (1|block), data = biomass, REML = F)

### check assumptions
check_model(r2)

#find optimal lambda for Box-Cox transformation 
#bc <- boxcox(ratio ~ treatment*variety, data = biomass)
#(lambda <- bc$x[which.max(bc$y)])
#fit new linear regression model using the Box-Cox transformation
#new_model <- lm(((ratio^lambda-1)/lambda) ~ treatment*variety, data = biomass)
#Q-Q plot for Box-Cox transformed model
qqnorm(new_model$residuals)
qqline(new_model$residuals)

AIC(r2, new_model)
AIC(r2) + 2*sum(log(biomass$ratio))
# r2         8 21.80522   4.2
# new_model  7 20.84490

Anova(r2)
#Anova(new_model)

br2 = emmeans(r2, ~ treatment*variety, type = "response", adjust = "none")
br2
#br = emmeans(new_model, ~ treatment*variety, type = "response", adjust = "none")
#br

cld(br2,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    #type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")    ### no adjustment = LSD for multiple comparisons


########### SLA (specific leaf area) #############
read.csv('SLA.csv')
sla <- read.csv('SLA.csv', header = T, sep = ",")

#make dates into a date format  
sla$date <- as.Date(sla$Date,'%m/%d/%Y')

#create DAT column 
sla$dat <- as.factor(difftime(as.POSIXct(sla$date), as.POSIXct('2021-06-14', tz="PT"), units="days"))

## remove dat 66 due to sample mislabeling!
sla <- subset(sla, dat != 66)

#factor variety
sla$variety <- as.factor(sla$Variety)

#factor treatment 
sla$treatment <- as.factor(sla$Treatment)

# create blocks for pot.number
# create blocks 
sla <- sla %>% 
  mutate(block = case_when(
    Row == "1" ~ "1",
    Row == "2" ~ "1", 
    Row == "3" ~ "2", 
    Row == "4" ~ "2", 
    Row == "5" ~ "3",
    Row == "6" ~ "3", 
    Row == "7" ~ "4", 
    Row == "8" ~ "4"))

###Check to see that all combinations are correct
xtabs(~ block + Plant.Number, data = sla)
xtabs(~ treatment + dat, data = sla)
xtabs(~ Plant.Number + dat, data = sla)
xtabs(~ treatment + variety, data = sla)

#Combined Variety & Treatment 
sla$vartreat <- as.factor(paste(sla$variety, sla$treatment))

plot(sla$SLA.cm.2.g ~ sla$vartreat)

library(ggplot2)
sla.stats.dat <- ddply(sla, c("variety", "treatment", 'dat'),
                       summarise,
                       N    = length(SLA.cm.2.g),#We use length instead of count.
                       mean = mean(SLA.cm.2.g),
                       sd   = sd(SLA.cm.2.g),
                       se   = sd / sqrt(N))
sla.stats.dat

line.sla.stats <- ggplot(sla.stats.dat, aes(x=dat, y=mean, group = treatment, color = treatment)) + 
  facet_wrap(~ variety) + 
  geom_point(size = 4) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, size=1,
                position=position_dodge(0)) +
  theme_classic() +
  ylab('SLA (cm^2/g)') +
  xlab('Days After Transplanting') +
  labs(color = 'Treatment', shape = 'Treatment') + #rename legend title
  theme(axis.title.x = element_text(size = 14), axis.text = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  #xlim('11', '18', '25', '32', '39') +
  ggtitle('SLA Field - 2021') +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title=element_text(size=14)) +
  theme(strip.text = element_text(size = 14))
line.sla.stats

library(plyr)
sla.stats <- ddply(sla, c("variety", "treatment"),
                   summarise,
                   N    = length(SLA.cm.2.g),#We use length instead of count.
                   mean = mean(SLA.cm.2.g),
                   sd   = sd(SLA.cm.2.g),
                   se   = sd / sqrt(N))
sla.stats

#Bar plot of mean and s.e.
sla.treatment.bar <- ggplot(sla.stats, aes(x = treatment, y=mean))+
  facet_wrap(~ variety) +
  #facet_grid(vars(Leaf.Position), vars(variety), scales = "free", space = 'free') + 
  geom_bar(stat = 'identity', position="dodge")+ 
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab("SLA (cm^2/g)") + 
  xlab("Treatment") +
  ggtitle('SLA - Field 2021') + 
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title=element_text(size=14)) + 
  theme(strip.text = element_text(size = 14)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,70)) + 
  #scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = expansion(add=c(0.2,0.2)))
#geom_text(aes(label=letters.chl, y = mean + se), vjust = -.5, size = 6)
sla.treatment.bar

library(lme4)
s1 <- lmer(SLA.cm.2.g ~ treatment*variety + (1|Row) + (1|dat), data = sla, REML = F)
s3 <- lmer(log(SLA.cm.2.g) ~ treatment*variety + (1|Row) + (1|dat), data = sla, REML = F)
s5 <- lmer(sqrt(SLA.cm.2.g) ~ treatment*variety + (1|Row) + (1|dat), data = sla, REML = F)

library(phia)
plot(interactionMeans(s1))
plot(interactionMeans(s3))
plot(interactionMeans(s5))

library(merTools)
RMSE.merMod(s1)
RMSE.merMod(s3)
RMSE.merMod(s5)

library(see)
library(performance)
### check assumptions
check_model(s1)
check_model(s3)
check_model(s5)

Anova(s1)
Anova(s3)
Anova(s5)

# estimated marginal means
ss2 = emmeans(s3, ~ treatment*variety, adjust = "none")
ss2

cld(ss2,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    #type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")    ### no adjustment = LSD for multiple comparisons

letters.sla = c("ab", "abc", "c", "bc", "a", "abc")

sla.stats <- ddply(sla, c("variety", "treatment"),
                   summarise,
                   N    = length(SLA.cm.2.g),#We use length instead of count.
                   mean = mean(SLA.cm.2.g),
                   sd   = sd(SLA.cm.2.g),
                   se   = sd / sqrt(N))
sla.stats

levels(sla.stats$variety)[levels(sla.stats$variety)=='Tahoe'] <- 'Tahoe Cinco'

library(ggplot2)
#Bar plot of mean and s.e.
sla.treatment.bar <- ggplot(sla.stats, aes(x = treatment, y=mean))+
  facet_wrap(~ variety) +
  #facet_grid(vars(Leaf.Position), vars(variety), scales = "free", space = 'free') + 
  geom_bar(stat = 'identity', position="dodge")+ 
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab("SLA (cm^2/g)") + 
  xlab("Treatment") +
  ggtitle('SLA - Field 2021') + 
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title=element_text(size=14)) + 
  theme(strip.text = element_text(size = 14)) + 
  scale_y_continuous(breaks=0:70*10, expand = c(0,0), limits = c(0,70)) + 
  scale_x_discrete(expand = expansion(add=c(0.2,0.2))) +
  geom_text(aes(label=letters.sla, y = mean + se), vjust = -.5, size = 6)
sla.treatment.bar


######### Leaf Chlorophyll Total #######
read.csv('Field Chlorophyll Data2.csv')
chlorophyll <- read.csv("Field Chlorophyll Data2.csv", header = T, sep = ",")

library(questionr)
#chlorophyll <- na.rm(chlorophyll)

#make dates into a date format  
chlorophyll$date <- as.Date(chlorophyll$Date,'%m/%d/%Y')

library(questionr)
#create DAT column 
chlorophyll$dat <- as.factor(difftime(as.POSIXct(chlorophyll$date), as.POSIXct('2021-06-14', tz="PT"), units="days"))

#factor variety
chlorophyll$variety <- as.factor(chlorophyll$Variety)

#factor treatment 
chlorophyll$treatment <- as.factor(chlorophyll$Treatment)

#factor block 
chlorophyll$block <- as.factor(chlorophyll$Block)

#plant number 
chlorophyll$plant.number <- as.factor(chlorophyll$Plant.Number)

# rename Tahoe Cinco
levels(chlorophyll$variety)[levels(chlorophyll$variety)=='Tahoe'] <- 'Tahoe Cinco'


#rename chlorophyll values
chlorophyll$chla.ug.ml <- chlorophyll$chla..ug.mL.
chlorophyll$chla.mg.g.fw <- chlorophyll$chla..mg.g.FW.

chlorophyll$chlb.ug.ml <- chlorophyll$chlb..ug.mL.
chlorophyll$chlb.mg.g.fw <- chlorophyll$chlb..mg.g.FW.

chlorophyll$chlt.ug.ml <- chlorophyll$chlt..ug.mL.
chlorophyll$chlt.mg.g.fw <- chlorophyll$chlt..mg.g..FW

# cm^2 chlorophyll per area
chlorophyll$chlacm2 <- chlorophyll$chla..mg.cm2.
chlorophyll$chlbcm2 <- chlorophyll$chlb.mg.cm2.
chlorophyll$chltcm2 <- chlorophyll$chlt..mg.cm2.

###Check to see that all combinations are correct
xtabs(~ variety + dat, data = chlorophyll)
xtabs(~ treatment + dat, data = chlorophyll)
xtabs(~ treatment + variety, data = chlorophyll)

library(rstatix)
select.chlorophyll <- chlorophyll %>% df_select(dat, plant.number, variety, treatment, block, chla.ug.ml, chla.mg.g.fw, chlb.ug.ml, chlb.mg.g.fw, 
                                                chlt.ug.ml, chlt.mg.g.fw, chlacm2, chlbcm2, chltcm2)

#Combined Variety & Treatment 
select.chlorophyll$vartreat <- as.factor(paste(select.chlorophyll$variety, select.chlorophyll$treatment))

library(plyr)
chl.stats <- ddply(select.chlorophyll, c("variety", "treatment"),
                   summarise,
                   N    = length(chltcm2),#We use length instead of count.
                   mean = mean(chltcm2),
                   sd   = sd(chltcm2),
                   se   = sd / sqrt(N))
chl.stats

ggplot(chl.stats, aes(x=treatment, y=mean, fill = treatment)) + 
  facet_wrap(~ variety) +
  geom_bar(stat = 'identity', position="dodge")+ 
  scale_fill_manual(values = c("#1B9E77", '#E7298A'))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, size=1,position=position_dodge(0)) +
  theme_classic() +
  ylab(bquote('Chlorophyll Total  '(mg~cm^-2))) +
  xlab('Treatment') +
  labs(color = 'Treatment', linetype = 'Cultivar') + #rename legend title
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title=element_text(size=14)) + 
  theme(strip.text = element_text(size = 14)) + 
  theme(axis.title.y = element_text(vjust=0.9)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.03)) + 
  scale_x_discrete(expand = expansion(add=c(0.2,0.2)))+
  #geom_text(aes(label=letters.chl, y = mean + se), vjust = -.5, size = 6)+
  theme(legend.position = 'none')

library(lme4)
### USE ML when comparing models with different fixed effects
lm1 <- lmer(chltcm2 ~ treatment*variety + (1|block) + (1|dat), data = select.chlorophyll, REML = F)
lm3 <- lmer(chltcm2 ~ treatment*variety + (1|block/plant.number) + (1|dat), data = select.chlorophyll, REML = F)
# model selection 
AIC(lm1, lm3)
#lm1  9 -1188.307
#lm3 10 -1189.531

#lm4 <- lmer(log(chltcm2) ~ treatment*variety + (1|block/plant.number) + (1|dat), data = select.chlorophyll, REML = F)
#lm5 <- lmer(sqrt(chltcm2) ~ treatment*variety + (1|block/plant.number) + (1|dat), data = select.chlorophyll, REML = F)

library(see)
library(performance)
### check assumptions
check_model(lm1)
check_model(lm3)

Anova(lm1)
#Anova(lm3)

library(multcompView)
library(emmeans)
library(multcomp)
# estimated marginal means
# e1 = emmeans(lm3, specs = pairwise ~ treatment|variety, adjust = "none", type = "response") # compare levels of treatment within group
# e1
e2 = emmeans(lm1, ~ treatment*variety, type = "response", adjust = "none")
e2
#e3 = emmeans(lm3, ~ treatment*variety, type = "response", adjust = "none")
#e3
# comparisons 
# confidence intervals with statistical tests 
contrast(e2, method="pairwise",adjust="none", type = "response", infer = T)

pairs(e2, alpha = 0.05, adjust = "none") # contrast stats without CI's
#write.csv(SPAD.contrasts, file = 'SPAD.contrasts.csv', quote = FALSE, row.names = T)

#summary(contrast(e2, method="pairwise",adjust="tukey"), infer=c(TRUE, TRUE)) # same as pairs

pwpp(e2, by = "treatment", adjust = "none")
pairs(e2, adjust = "none")
plot(e2, comparisons = T, adjust = "none", by = "variety") # if red arrows overlap groups are not sig from each other

#contrast(regrid(e2))
cld(e2,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", ### no adjustment = LSD for multiple comparisons
    method = "paiwise")    

letters.chl = c("c", "b", "c", "ab", "b", "a")

### manuscript chlorphyll figure 
chl.treatment.bar.bw <- ggplot(chl.stats, aes(x = treatment, y=mean, fill = treatment))+
  facet_wrap(~ variety) +
  #facet_grid(vars(Leaf.Position), vars(variety), scales = "free", space = 'free') + 
  geom_bar(stat = 'identity', position="dodge", color = "black")+ 
  scale_fill_grey()+
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab(bquote('Chlorophyll Total  '(mg~cm^-2))) +
  xlab("N Treatment") +
  #ggtitle('SPAD and Chlorophyll - Field 2021') + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  theme(legend.text=element_text(size=16)) + 
  theme(legend.title=element_text(size=16)) + 
  theme(strip.text = element_text(size = 16)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,0.03)) + 
  theme(legend.position = 'none') +
  geom_text(aes(label=letters.chl, y = mean + se), vjust = -.5, size = 7) 
chl.treatment.bar.bw


################## Leaf % Nitrogen and delta Carbon 13 Data ###############################
read.csv('nitrogen and carbon data.csv')
nc <- read.csv("nitrogen and carbon data.csv", header = T, sep = ",")

library(questionr)
#make dates into a date format  
nc$date <- as.Date(nc$Date,'%m/%d/%Y')

#create DAT column 
nc$dat <- as.factor(difftime(as.POSIXct(nc$date), as.POSIXct('2021-06-14', tz="PT"), units="days"))

#factor variety
nc$variety <- as.factor(nc$Variety)

# rename Tahoe Cinco
levels(nc$variety)[levels(nc$variety)=='Tahoe'] <- 'Tahoe Cinco'

#factor treatment 
nc$treatment <- as.factor(nc$Treatment)

#factor block 
nc$block <- as.factor(nc$Block)

#plant number 
nc$plant.number <- as.factor(nc$Plant.Number)

#Combined Variety & Treatment 
nc$vartreat <- as.factor(paste(nc$variety, nc$treatment))

##### Leaf Nitrogen % ########
n.stats <- ddply(nc, c("variety", "treatment"),
                 summarise,
                 N    = length(perc.Nitrogen),#We use length instead of count.
                 mean = mean(perc.Nitrogen),
                 sd   = sd(perc.Nitrogen),
                 se   = sd / sqrt(N))
n.stats

#Bar plot of mean and s.e.
n.bar <- ggplot(n.stats, aes(x = treatment, y=mean))+
  facet_wrap(~ variety) +
  #facet_grid(vars(Leaf.Position), vars(variety), scales = "free", space = 'free') + 
  geom_bar(stat = 'identity', position="dodge")+ 
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab("Leaf Nitrogen %") + 
  xlab("Treatment") +
  #ggtitle('Phi NPQ - Field 2021') + 
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title=element_text(size=14)) + 
  theme(strip.text = element_text(size = 11)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,5.5))
  #geom_text(aes(label=letters.percn, y = mean + se), vjust = -.5, size = 6)
n.bar

m1 <- lmer(perc.Nitrogen ~ treatment*variety + (1|block) + (1|dat), data = nc, REML = F)
m4 <- lmer(perc.Nitrogen ~ treatment+variety + (1|block) + (1|dat), data = nc, REML = F)

m7 <- lmer(sqrt(perc.Nitrogen) ~ treatment*variety + (1|block) + (1|dat), data = nc, REML = F)
m10 <- lmer(sqrt(perc.Nitrogen) ~ treatment+variety + (1|block) + (1|dat), data = nc, REML = F)

m13 <- lmer(log(perc.Nitrogen) ~ treatment*variety + (1|block) + (1|dat), data = nc, REML = F)
m16 <- lmer(log(perc.Nitrogen) ~ treatment+variety + (1|block) + (1|dat), data = nc, REML = F)

library(phia)
plot(interactionMeans(m7))
plot(interactionMeans(m10))
plot(interactionMeans(m13))
plot(interactionMeans(m16))

library(merTools)
RMSE.merMod(m1)
RMSE.merMod(m4)
RMSE.merMod(m7)
RMSE.merMod(m10)
RMSE.merMod(m13)
RMSE.merMod(m16)

library(see)
library(performance)
### check assumptions
check_model(m7)
check_model(m10)
check_model(m13)
check_model(m16)

#Anova(m16)
Anova(m13)

# estimated marginal means
e2 = emmeans(m13, ~ treatment*variety, adjust = "none", type = "response")
e2

cld(e2,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    #type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", ### no adjustment = LSD for multiple comparisons
    method = "paiwise") 

letters.percn = c("cd", "ab", "d", "ab", "bc", "a")

percn.var <- emmeans(e2, specs = pairwise ~ variety)

cld(percn.var, alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    #type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")

### Leaf Nitrogen % manuscript figure  
n.bar.bw <- ggplot(n.stats, aes(x = treatment, y=mean, fill = treatment))+
  facet_wrap(~ variety) +
  #facet_grid(vars(Leaf.Position), vars(variety), scales = "free", space = 'free') + 
  geom_bar(stat = 'identity', position="dodge", color = "black")+ 
  scale_fill_grey()+
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab("Leaf Nitrogen (%)") + 
  xlab("Treatment") +
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 18)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  theme(legend.text=element_text(size=18)) + 
  theme(legend.title=element_text(size=18)) + 
  theme(strip.text = element_text(size = 18)) + 
  #scale_y_continuous(breaks=0:70*10, expand = c(0,0), limits = c(0,70)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 6)) +
  theme(legend.position = 'none') +
  geom_text(aes(label=letters.percn, y = mean + se), vjust = -.5, size = 7) 
n.bar.bw

######### delta 13C #######
deltac.stats <- ddply(nc, c("variety", "treatment"),
                      summarise,
                      N    = length(delta.13.C),#We use length instead of count.
                      mean = mean(delta.13.C),
                      sd   = sd(delta.13.C),
                      se   = sd / sqrt(N))
deltac.stats

#Bar plot of mean and s.e.
deltac.bar <- ggplot(deltac.stats, aes(x = treatment, y=mean))+
  facet_wrap(~ variety) +
  #facet_grid(vars(Leaf.Position), vars(variety), scales = "free", space = 'free') + 
  geom_bar(stat = 'identity', position="dodge")+ 
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab("delta 13 C") + 
  xlab("Treatment") +
  #ggtitle('Phi NPQ - Field 2021') + 
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title=element_text(size=14)) + 
  theme(strip.text = element_text(size = 14)) +
  scale_y_continuous(expand = c(0,0), limits = c(-30, 0))
  #geom_text(aes(label=letters.deltac, y = mean + se), vjust = 1.9, size = 6)
deltac.bar

interaction.plot(x.factor     = nc$treatment,
                 trace.factor = nc$variety,
                 response     = nc$delta.13.C,
                 fun = mean,
                 type="b",
                 col=c("black","red","green"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

dc1 <- lmer(delta.13.C ~ treatment*variety + (1|block) + (1|dat), data = nc, REML = F)
dc4 <- lmer(delta.13.C ~ treatment+variety + (1|block) + (1|dat), data = nc, REML = F)
#dc5 <- lmer(delta.13.C ~ treatment+variety + (1|block/plant.number) + (1/dat), data = nc, REML = F)
dclog <- lmer(log(-1* delta.13.C) ~ treatment*variety + (1|block) + (1|dat), data = nc, REML = F)

AIC(dc1, dc4) 
# df      AIC
# dc1  9 386.1124
#dc4  7 398.2950

library(see)
library(performance)
### check assumptions
check_model(dc1)
check_model(dclog)

library(merTools)
RMSE.merMod(dc1)
RMSE.merMod(dclog)

Anova(dclog)
Anova(dc1)

library(phia)
plot(interactionMeans(dc1))

library(multcompView)
library(emmeans)
library(multcomp)
# estimated marginal means
e1 = emmeans(dc1, ~ treatment*variety, adjust = "none", type = "response") 
e1

cld(e1,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", ### no adjustment = LSD for multiple comparisons
    method = "paiwise") 

letters.deltac = c("d", "bc", "c", "a", "ab", "a")

### black and white 
deltac.bar.bw <- ggplot(deltac.stats, aes(x = treatment, y=mean, fill = treatment))+
  facet_wrap(~ variety) +
  #facet_grid(vars(Leaf.Position), vars(variety), scales = "free", space = 'free') + 
  geom_bar(stat = 'identity', position="dodge", color = "black")+ 
  scale_fill_grey()+
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab(bquote(~delta^13*'C')) +
  #ylab("Phi NPQ") + 
  xlab("N Treatment") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  theme(legend.text=element_text(size=16)) + 
  theme(legend.title=element_text(size=16)) + 
  theme(strip.text = element_text(size = 16)) +
  coord_cartesian(ylim=c(-30,-20)) +
  scale_y_continuous(breaks=-30:0*2, expand = c(0,0), limits = c(-30,0)) + 
  #scale_y_continuous(expand = c(0,0))+
  theme(legend.position = 'none') +
  geom_text(aes(label=letters.deltac, y = mean + se), vjust = 2, size = 7, 
            nudge_y = c(-0.4, -0.1, -0.1, -0.1, -0.1, 0)) 
deltac.bar.bw

#### merge Leaf Nitrogen % & delta 13C manuscript figures 
leafn.deltac.bw2 <- plot_grid(n.bar.bw, deltac.bar.bw, labels = "AUTO", label_size = 18)
leafn.deltac.bw2
#ggsave(leafn.deltac.bw2, file="leafn.deltac.bw.TIFF",width=12, height=5,dpi=600,path="C:/Users/mfarnisa/OneDrive - University of Nevada, Reno/Thesis/Figures")



########### SLN (specific leaf nitrogen) ########
short.sla <- sla %>% df_select(dat, Plant.Number, variety, treatment, block, SLA.cm.2.g)
short.nc <- nc %>% df_select(dat, Plant.Number, variety, treatment, block, perc.Nitrogen)

slnjoin <- full_join(short.nc, short.sla, by = c('dat', 'variety', 'treatment', 'block', 'Plant.Number'))

## remove dat 66 due to sample mislabeling!
slnjoin <- subset(slnjoin, dat != 66)

# calculate N ug/cm2
slnjoin$SLNug.cm2 = ((1/slnjoin$SLA.cm.2.g)*(slnjoin$perc.Nitrogen/100) * 1000000)

sln.stats <- ddply(slnjoin, c("variety", "treatment"),
                   summarise,
                   N    = length(SLNug.cm2),#We use length instead of count.
                   mean = mean(SLNug.cm2),
                   sd   = sd(SLNug.cm2),
                   se   = sd / sqrt(N))
sln.stats

library(lme4)
s1 <- lmer(SLNug.cm2 ~ treatment*variety + (1|block) + (1|dat), data = slnjoin, REML = F)
s1.5 <- lmer(SLNug.cm2 ~ treatment+variety + (1|block) + (1|dat), data = slnjoin, REML = F)

s3 <- lmer(log(SLNug.cm2) ~ treatment*variety + (1|block) + (1|dat), data = slnjoin, REML = F)
s4 <- lmer(log(SLNug.cm2) ~ treatment+variety + (1|block) + (1|dat), data = slnjoin, REML = F)

s5 <- lmer(sqrt(SLNug.cm2) ~ treatment*variety + (1|block) + (1|dat), data = slnjoin, REML = F)
s6 <- lmer(sqrt(SLNug.cm2) ~ treatment+variety + (1|block) + (1|dat), data = slnjoin, REML = F)

library(phia)
plot(interactionMeans(s3))
plot(interactionMeans(s4))

library(merTools)
RMSE.merMod(s3)
RMSE.merMod(s4)
RMSE.merMod(s5)
RMSE.merMod(s6)

library(see)
library(performance)
### check assumptions
check_model(s3)
check_model(s4)

Anova(s3)
Anova(s4)

library(multcompView)
library(emmeans)
library(multcomp)
# estimated marginal means

ss3 = emmeans(s3, ~ treatment*variety, adjust = "none")
ss3

cld(ss3,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    #type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")    ### no adjustment = LSD for multiple comparisons

letters.sln = c("b", "ab", "ab", "a", "b", "ab")

#Bar plot of mean and s.e.
sln.bar <- ggplot(sln.stats, aes(x = treatment, y=mean))+
  facet_wrap(~ variety) +
  #facet_grid(vars(Leaf.Position), vars(variety), scales = "free", space = 'free') + 
  geom_bar(stat = 'identity', position="dodge")+ 
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab(bquote('SLN  '(ug~cm^-2))) +
  #ylab("SLN (ug/cm2)") + 
  xlab("Treatment") +
  #ggtitle('Phi NPQ - Field 2021') + 
  theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 16)) +
  theme(legend.text=element_text(size=16)) + 
  theme(legend.title=element_text(size=16)) + 
  theme(strip.text = element_text(size = 16)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1200))+
  geom_text(aes(label=letters.sln, y = mean + se), vjust = -.5, size = 6)
sln.bar


############### SPAD - PhotosynQ Measurements ############ 
read.csv('photosynq data.csv')
pq <- read.csv("photosynq data.csv", header = T, sep = ",")

#factor variety
pq$variety <- as.factor(pq$Variety)
pq$treatment <- as.factor(pq$Treatment)
pq$plant.number <- as.factor(pq$Plant.Number)
pq$row <- as.factor(pq$Row)

# rename Tahoe Cinco
levels(pq$variety)[levels(pq$variety)=='Tahoe'] <- 'Tahoe Cinco'

#make dates into a date format  
pq$date <- as.Date(pq$Date,'%m/%d/%Y')

library(questionr)
#create DAT column 
pq$dat <- as.factor(difftime(as.POSIXct(pq$date), as.POSIXct('2021-06-14', tz="PT"), units="days"))

allpq <- pq[c('SPAD', "FvP_over_FmP","Phi2", "variety", "treatment", "plant.number", "row", "dat", "date")]

# create blocks 
allpq <- allpq %>% 
  mutate(block = case_when(
    row == "1" ~ "1",
    row == "2" ~ "1", 
    row == "3" ~ "2", 
    row == "4" ~ "2", 
    row == "5" ~ "3",
    row == "6" ~ "3", 
    row == "7" ~ "4", 
    row == "8" ~ "4"))

allpq$block <- as.factor(allpq$block)
allpq$vartreat <- as.factor(paste(allpq$variety, allpq$treatment))

#calculate daily average for photosynq measurements  
pqavgs1 <- aggregate(cbind(FvP_over_FmP, SPAD, Phi2) ~ dat + plant.number + variety + treatment + vartreat + block, data=allpq, mean)

xtabs(~ variety + dat, data = spadavgs)
xtabs(~ variety + dat, data = pqavgs1)
xtabs(~ treatment + dat, data = spadavgs)
xtabs(~ plant.number + dat, data = spadavgs)
xtabs(~ treatment + variety, data = spadavgs)

plot(mean ~ dat, data = spadavgs)
plot(SPAD ~ dat, data = pqavgs1)

m1 <- lmer(SPAD ~ treatment*variety + (1|block) + (1|dat), data = pqavgs1, REML = F)
summary(m1)
# check plots 
plot(m1)
par(mfrow=c(1,2))
qqnorm(resid(m1))
qqline(resid(m1))
# residuals
hist(residuals(m1), breaks=10)
anova(m1)


m2 <- lmer(SPAD ~ treatment*variety + (1|block/plant.number) + (1/dat), data = pqavgs1, REML = F)
summary(m2)
# check plots
plot(m2)
par(mfrow=c(1,2))
qqnorm(resid(m2))
qqline(resid(m2))
# residuals
hist(residuals(m2), breaks=10)
anova(m2)

m3 <- lmer(SPAD ~ treatment*variety + (1|plant.number) + (1|dat), data = pqavgs1, REML = F)
summary(m3)
# check plots 
plot(m3)
#par(mfrow=c(1,2))
qqnorm(resid(m3))
qqline(resid(m3))
# residuals
hist(residuals(m3), breaks=10)
anova(m3)

m4 <- lmer(SPAD ~ treatment+variety + (1|block) + (1|dat), data = pqavgs1, REML = F)
summary(m4)
# check plots
plot(m4)
par(mfrow=c(1,2))
qqnorm(resid(m4))
qqline(resid(m4))
# residuals
hist(residuals(m4), breaks=10)

AIC(m1, m2, m3, m4) #m3
# df      AIC
# m1  9 6107.921
# m2  9 6324.435
# m3  9 6092.669
# m4  7 6123.806

anova(m3, m1, m4)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# m4    7 6123.8 6155.9 -3054.9   6109.8                         
# m3    9 6092.7 6133.9 -3037.3   6074.7 35.137  2  2.345e-08 ***
# m1    9 6107.9 6149.1 -3045.0   6089.9  0.000  0               

library(see)
library(performance)
### check assumptions
check_model(m3)

# estimated marginal means
e1 = emmeans(m3, specs = pairwise ~ treatment|variety, adjust = "none", type = "response") # compare levels of treatment within group
e1
e2 = emmeans(m3, ~ treatment*variety, adjust = "none", type = "response")
e2
e3 = emmeans(m3,  ~ treatment|variety, type = "response") # compare levels of treatment within group
e3
# comparisons 
# confidence intervals with statistical tests 
contrast(e3, method="pairwise",adjust="none", type = "response", infer = T)

SPAD.contrasts <- pairs(e3, alpha = 0.05) # contrast stats without CI's
#write.csv(SPAD.contrasts, file = 'SPAD.contrasts.csv', quote = FALSE, row.names = T)
#summary(contrast(e2, method="pairwise",adjust="tukey"), infer=c(TRUE, TRUE)) # same as pairs

pwpp(e2, by = "treatment")

plot(e2, comparisons = T) # if red arrows overlap groups are not sig from each other

cld(e2)
SPAD.lmer.cld <- cld(e2,
                     alpha   = 0.05,
                     reversed  = T, ### reverse the order of letters
                     Letters = letters,    ### Use lower-case letters for .group
                     type    = "response", ### Report emmeans in orginal scale
                     adjust =  "none", ### no adjustment = LSD for multiple comparisons
                     method = "paiwise") 
SPAD.lmer.cld
#write.csv(SPAD.lmer.cld, file = 'SPAD.cld.csv', quote = FALSE, row.names = T)

letters.spad = c("c", "b", "c", "ab", "a", "a")

library(plyr)
spad.vis <- ddply(pqavgs1, c("variety", "treatment"),
                  summarise,
                  N    = length(SPAD),#We use length instead of count.
                  mean = mean(SPAD),
                  sd   = sd(SPAD),
                  se   = sd / sqrt(N))
spad.vis

### manuscript SPAD figure 
spad.treatment.bar.bw <- ggplot(spad.vis, aes(x = treatment, y=mean, fill = treatment))+
  facet_wrap(~ variety) +
  #facet_grid(vars(Leaf.Position), vars(variety), scales = "free", space = 'free') + 
  geom_bar(stat = 'identity', position="dodge", color = "black")+ 
  scale_fill_grey()+
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab("SPAD") + 
  xlab("Treatment") +
  #ggtitle('SPAD and Chlorophyll - Field 2021') + 
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 18)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  theme(legend.text=element_text(size=18)) + 
  theme(legend.title=element_text(size=18)) + 
  theme(strip.text = element_text(size = 18)) + 
  scale_y_continuous(breaks=0:75*10, expand = c(0,0), limits = c(0,75)) + 
  theme(legend.position = 'none') +
  geom_text(aes(label=letters.spad, y = mean + se), vjust = -.5, size = 7) 
spad.treatment.bar.bw
#ggsave(spad.treatment.bar.bw, file="field22 spad.bw.TIFF",width=13, height=6,dpi=600,path="H:/WNM Conference 2023")


########### Fv'/Fm' - PhotosynQ ############
fvfm.stats.dat <- ddply(pqavgs1, c("vartreat", "dat"),
                        summarise,
                        N    = length(FvP_over_FmP),#We use length instead of count.
                        mean = mean(FvP_over_FmP),
                        sd   = sd(FvP_over_FmP),
                        se   = sd / sqrt(N))
fvfm.stats.dat

line.fvfm.stats <- ggplot(fvfm.stats.dat, aes(x=dat, y=mean, group = vartreat, color = vartreat, shape = vartreat)) + 
  geom_point(size = 4) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, size=1,
                position=position_dodge(0)) +
  theme_classic() +
  ylab('FvP/FmP') +
  xlab('Days After Transplanting') +
  labs(color = 'Treatment', shape = 'Treatment') + #rename legend title
  theme(axis.title.x = element_text(size = 14), axis.text = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  #scale_y_continuous(breaks=0:1*.2) +
  ggtitle('PhotosynQ Field 2021') +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title=element_text(size=14)) +
  theme(strip.text = element_text(size = 14))
line.fvfm.stats

# remove dat 80
fvfm.reduced <- subset(pqavgs1, dat != 64 & dat != 65 &dat != 66 &dat != 70 & dat != 80 &dat != 81& dat != 87 &dat != 88)
fvfm.reduced1 <- subset(pqavgs1, dat != 65 & dat != 66  & dat != 64  &dat != 70)

# check outliers
(pqavgs1$FvP_over_FmP < 0) == T
sort(pqavgs1$FvP_over_FmP)

library(plyr)
fvfm.vis <- ddply(pqavgs1, c("variety", "treatment"),
                  summarise,
                  N    = length(FvP_over_FmP),#We use length instead of count.
                  mean = mean(FvP_over_FmP),
                  sd   = sd(FvP_over_FmP),
                  se   = sd / sqrt(N))
fvfm.vis

library(ggplot2)
#Bar plot of mean and s.e.
fvfm.treatment.bar <- ggplot(fvfm.vis, aes(x = treatment, y=mean))+
  facet_wrap(~ variety) +
  #facet_grid(vars(Leaf.Position), vars(variety), scales = "free", space = 'free') + 
  geom_bar(stat = 'identity', position="dodge")+ 
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab("FvP/FmP") + 
  xlab("Treatment") +
  #ggtitle('PhotosynQ Field 2021') + 
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  #theme(legend.text=element_text(size=12)) + 
  #theme(legend.title=element_text(size=14)) + 
  theme(strip.text = element_text(size = 14)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.8))
  #scale_x_discrete(expand = expansion(add=c(0.2,0.2))) 
fvfm.treatment.bar

interaction.plot(x.factor     = pqavgs1$treatment,
                 trace.factor = pqavgs1$variety,
                 response     = pqavgs1$FvP_over_FmP,
                 fun = mean,
                 type="b",
                 col=c("black","red","green"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

m1 <- lmer(FvP_over_FmP ~ treatment*variety + (1|block/plant.number) + (1|dat), data = pqavgs1, REML = F)
m4 <- lmer(FvP_over_FmP ~ treatment+variety + (1|block/plant.number) + (1|dat), data = pqavgs1, REML = F)

#m7 <- lmer(sqrt(FvP_over_FmP) ~ treatment*variety + (1|block/plant.number) + (1|dat), data = pqavgs1, REML = F)
#m10 <- lmer(sqrt(FvP_over_FmP) ~ treatment+variety + (1|block/plant.number) + (1|dat), data = pqavgs1, REML = F)
#m11 <- lmer(log(FvP_over_FmP) ~ treatment*variety + (1|block/plant.number) + (1|dat), data = pqavgs1, REML = F)
#m12 <- lmer(log(FvP_over_FmP) ~ treatment+variety + (1|block/plant.number) + (1|dat), data = pqavgs1, REML = F)

library(rcompanion)
T_tuk = transformTukey(pqavgs1$FvP_over_FmP, plotit=T)

fvfm.sqr = (pqavgs1$FvP_over_FmP)^(4.975)   # Avoid complex numbers

m13 <- lmer(fvfm.sqr ~ treatment*variety + (1|block/plant.number) + (1|dat), data = pqavgs1, REML = F)
m16 <- lmer(fvfm.sqr ~ treatment+variety + (1|block/plant.number) + (1|dat), data = pqavgs1, REML = F)


library(see)
library(performance)
### check assumptions
check_model(m1)  
check_model(m4)  
# check_model(m7)  
# check_model(m10)  
# check_model(m11)  
# check_model(m12)  
check_model(m13)
check_model(m16)  

library(merTools)
RMSE.merMod(m1)
RMSE.merMod(m4)
# RMSE.merMod(m7)
# RMSE.merMod(m10)
# RMSE.merMod(m11)
# RMSE.merMod(m12)
RMSE.merMod(m13)
RMSE.merMod(m16)

Anova(m1)
Anova(m4)
Anova(m13)
Anova(m16)

library(multcompView)
library(emmeans)
library(multcomp)
# estimated marginal means
mm5 = emmeans(m13, ~ treatment*variety, adjust = "none", type = "response") 
mm5
#mm2 = emmeans(m4, ~ treatment*variety, adjust = "none", type = "response")
#mm2

cld(mm5,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    #type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", ### no adjustment = LSD for multiple comparisons
    method = "paiwise") 


letters.fvfm = c("b", "ab", "ab", "ab", "a", "a")

### Fv'/Fm' manuscript figure 
fvfm.treatment.bar.bw <- ggplot(fvfm.vis, aes(x = treatment, y=mean, fill = treatment))+
  facet_wrap(~ variety) +
  #facet_grid(vars(Leaf.Position), vars(variety), scales = "free", space = 'free') + 
  geom_bar(stat = 'identity', position="dodge", color = "black")+ 
  scale_fill_grey()+
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab("Fv`/Fm`") + 
  xlab("N Treatment") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  theme(legend.text=element_text(size=16)) + 
  theme(legend.title=element_text(size=16)) + 
  theme(strip.text = element_text(size = 16)) + 
  #scale_y_continuous(breaks=0:70*10, expand = c(0,0), limits = c(0,70)) + 
  #ylim(0.5,0.8)+
  coord_cartesian(ylim=c(0.5,0.8)) +
  #scale_y_continuous(expand = c(0.5,0), limits = c(0.5, 0.8)) +
  theme(legend.position = 'none') +
  geom_text(aes(label=letters.fvfm, y = mean + se), vjust = -.5, size = 7) 
fvfm.treatment.bar.bw



############ Phi2 - PhotosynQ #############
plot(Phi2 ~ vartreat, data = pqavgs1)
sort(pqavgs1$Phi2)

library(plyr)
phi2.vis <- ddply(pqavgs1, c("variety", "treatment"),
                  summarise,
                  N    = length(Phi2),#We use length instead of count.
                  mean = mean(Phi2),
                  sd   = sd(Phi2),
                  se   = sd / sqrt(N))
phi2.vis
#write.table(phi2.vis, file = 'Phi2 field 2021 Means.csv', sep = ",", quote = FALSE, row.names = F)

ggplot(phi2.vis, aes(x=treatment, y=mean, fill = treatment)) + 
  facet_wrap(~ variety) +
  geom_bar(stat = 'identity', position="dodge")+ 
  scale_fill_manual(values = c("#1B9E77", '#E7298A'))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, size=1,position=position_dodge(0)) +
  theme_classic() +
  ylab(bquote(~phi*'PSII')) +
  xlab('Treatment') +
  labs(color = 'Treatment', linetype = 'Cultivar') + #rename legend title
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title=element_text(size=14)) + 
  theme(strip.text = element_text(size = 14)) + 
  theme(axis.title.y = element_text(vjust=0.9)) +
  coord_cartesian(ylim=c(0.3,0.6)) +
  scale_x_discrete(expand = expansion(add=c(0.2,0.2)))+
  #geom_text(aes(label=letters.phi2, y = mean + se), vjust = -.5, size = 6)+
  theme(legend.position = 'none')

library(ggplot2)
#Bar plot of mean and s.e.
phi2.treatment.bar <- ggplot(phi2.vis, aes(x = treatment, y=mean))+
  facet_wrap(~ variety) +
  #facet_grid(vars(Leaf.Position), vars(variety), scales = "free", space = 'free') + 
  geom_bar(stat = 'identity', position="dodge")+ 
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab("Phi2") + 
  xlab("Treatment") +
  #ggtitle('PhotosynQ Field 2021') + 
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  #theme(legend.text=element_text(size=12)) + 
  #theme(legend.title=element_text(size=14)) + 
  theme(strip.text = element_text(size = 14)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.6))
  #scale_x_discrete(expand = expansion(add=c(0.2,0.2))) +
  #geom_text(aes(label=letters.phi2, y = mean + se), vjust = -.5, size = 6)
phi2.treatment.bar

interaction.plot(x.factor     = pqavgs1$treatment,
                 trace.factor = pqavgs1$variety,
                 response     = pqavgs1$Phi2,
                 fun = mean,
                 type="b",
                 col=c("black","red","green"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

m1 <- lmer(Phi2 ~ treatment*variety + (1|block) + (1|dat), data = pqavgs1, REML = F)
m4 <- lmer(Phi2 ~ treatment+variety + (1|block) + (1|dat), data = pqavgs1, REML = F)

m7 <- lmer(sqrt(Phi2) ~ treatment*variety + (1|block) + (1|dat), data = pqavgs1, REML = F)
m10 <- lmer(sqrt(Phi2) ~ treatment+variety + (1|block) + (1|dat), data = pqavgs1, REML = F)

m13 <- lmer(log(Phi2) ~ treatment*variety + (1|block) + (1|dat), data = pqavgs1, REML = F)
m16 <- lmer(log(Phi2) ~ treatment+variety + (1|block) + (1|dat), data = pqavgs1, REML = F)

m17 <- lmer(phi2.sqr ~ treatment*variety + (1|block) + (1|dat), data = pqavgs1, REML = F)
m18 <- lmer(phi2.sqr ~ treatment+variety + (1|block) + (1|dat), data = pqavgs1, REML = F)

library(rcompanion)
phi2.transfrom = transformTukey(pqavgs1$Phi2, plotit=T)

phi2.sqr = (pqavgs1$Phi2)^(1.475)   # Avoid complex numbers

library(merTools)
RMSE.merMod(m1)
RMSE.merMod(m4)
RMSE.merMod(m7)
RMSE.merMod(m10)
RMSE.merMod(m13)
RMSE.merMod(m16)
RMSE.merMod(m17)
RMSE.merMod(m18)

library(see)
library(performance)
### check assumptions
check_model(m1)
check_model(m4)  
check_model(m7)  
check_model(m10)  
#check_model(m13)  
#check_model(m16)
check_model(m17)
check_model(m18)

Anova(m1)
Anova(m4)
Anova(m7)
Anova(m10)
Anova(m17)
Anova(m18)

# estimated marginal means
mm17 = emmeans(m17, ~ treatment*variety, adjust = "none", type = "response") 
mm17


cld(mm17,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    #type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", ### no adjustment = LSD for multiple comparisons
    method = "paiwise") 


letters.phi2 = c("bc", "abc", "a", "c", "ab", "abc")

### Phi2 manuscript figure 
phi2.treatment.bar.bw <- ggplot(phi2.vis, aes(x = treatment, y=mean, fill = treatment))+
  facet_wrap(~ variety) +
  #facet_grid(vars(Leaf.Position), vars(variety), scales = "free", space = 'free') + 
  geom_bar(stat = 'identity', position="dodge", color = "black")+ 
  scale_fill_grey()+
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab(bquote(~phi*'PSII')) +
  xlab("N Treatment") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  theme(legend.text=element_text(size=16)) + 
  theme(legend.title=element_text(size=16)) + 
  theme(strip.text = element_text(size = 16)) + 
  #scale_y_continuous(breaks=0:70*10, expand = c(0,0), limits = c(0,70)) + 
  #ylim(0.5,0.8)+
  coord_cartesian(ylim=c(0.3,0.6)) +
  #scale_y_continuous(expand = c(0.5,0), limits = c(0.5, 0.8)) +
  theme(legend.position = 'none') +
  geom_text(aes(label=letters.phi2, y = mean + se), vjust = -.5, size = 7) 
#labs(tag = "B") +
#theme(plot.tag = element_text(), plot.tag.position = c(0, 1))
phi2.treatment.bar.bw



######### Correlation Matrix #########
# merge data frames 
select.nc <- nc %>% df_select(dat, plant.number, variety, treatment, block, perc.Nitrogen, delta.13.C)
test <- left_join(pqavgs1, select.chlorophyll, by = c('dat', 'variety', 'treatment', 'block', 'plant.number'))
joined <- full_join(test, select.nc, by = c('dat', 'variety', 'treatment', 'block', 'plant.number'))

#Combined Variety & Treatment 
joined$vartreat <- as.factor(paste(joined$variety, joined$treatment))

# remove na's
joined1 <- na.rm(joined)

#factor plant.number
short.sla$plant.number <- short.sla$Plant.Number
short.sla$plant.number <- as.factor(short.sla$plant.number)

# merge data frames 
joined2 <- full_join(joined1, short.sla, by = c('dat', 'variety', 'treatment', 'block', 'plant.number'))

# calculate N ug/cm2
joined2$SLNug.cm2 = ((1/joined2$SLA.cm.2.g)*(joined2$perc.Nitrogen/100) * 1000000)

# improved correlation matrix
joined.corr <- joined2 %>% df_select(dat, variety, treatment, SPAD, FvP_over_FmP, Phi2, chltcm2, perc.Nitrogen, SLA.cm.2.g, SLNug.cm2)

colnames(joined.corr) <- c('dat', "variety", "treatment" , "SPAD", "Fv/Fm", "phi PSII", "total leaf chlorophyll", "total leaf N", "SLA", "SLN")

# change column names 
colnames(joined.corr)[5] <- 'Fv`/Fm`'
colnames(joined.corr)[6] <- 'Phi PSII'
colnames(joined.corr)[7] <- 'Total leaf chlorophyll'
colnames(joined.corr)[8] <- 'Total leaf N'

library(Hmisc)
library(corrplot)

parameters_rcorr <-as.matrix(joined.corr[,4:10])
paramaters_mat <-rcorr(parameters_rcorr)
# mat_2 <-rcorr(as.matrix(data)) returns the same output
paramaters_mat

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

cor.matrix <- corrplot(paramaters_mat$r, method="circle", col=col(200), 
                       type="lower", order="hclust", bg="white", addgrid.col = T,
                       addCoef.col = "black", # Add coefficient of correlation
                       tl.col="black", tl.srt=50, #Text label color and rotation
                       # Combine with significance
                       p.mat = paramaters_mat$P, sig.level = 0.05, insig = "blank", # hide correlation coefficient on the principal diagonal
                       diag=FALSE)

#ggsave(cor.matrix, file="cor.matrix.TIFF",width=12, height=5,dpi=600,path="C:/Users/mfarnisa/OneDrive - University of Nevada, Reno/Thesis/Figures")


############ Plants Flowering ###########
library(glmer)
library(glmmTMB)
flowering <- read.csv('comb.plantsflowering.csv')

flowering$variety <- as.factor(flowering$Variety)
flowering$treatment <- as.factor(flowering$Treatment)
flowering$vartreat <- as.factor(paste(flowering$variety, flowering$treatment))
flowering$dat <- as.factor(flowering$DAT)

levels(flowering$variety)[levels(flowering$variety)=='Tahoe '] <- 'Tahoe Cinco'

###### set up data frame like this example 
# Rootexod05.perc <- ddply(Rootanatomy, c("Rootstock","Treatment"),
#                          summarise,
#                          N    = length(Exod05),#We use length instead of count.
#                          sum = sum(Exod05),
#                          perc   = sum/N0)
# Rootexod05.perc
#modExo05<-glm(cbind(sum, N-sum)~Rootstock*Treatment,family=binomial,data=Rootexod05.perc)

m1 <- glm(cbind(Plants.Flowering,Not.flowering) ~ treatment, family = binomial,data = flowering)
m2 <- glm(cbind(Plants.Flowering,Not.flowering) ~ treatment*variety, family = binomial,data = flowering)
m3 <- glm(cbind(Plants.Flowering,Not.flowering) ~ treatment+variety, family = binomial,data = flowering)

m4 <- glmer(cbind(Plants.Flowering,Not.flowering) ~ treatment*variety + (1|dat), family = binomial,data = flowering)
#m5 <- glmer(cbind(Plants.Flowering,Not.flowering) ~ treatment+variety + (1|dat), family = binomial,data = flowering)
#m6 <- glm(cbind(Plants.Flowering,Not.flowering) ~ treatment*variety+dat, family = binomial,data = flowering)

AIC(m1, m2, m3)

Anova(m1)
Anova(m2)
Anova(m3)

AIC(m1, m2, m3, m4)

Anova(m4)

(m2.1 <- deviance(m2))
(m4.1 <- deviance(m4))

library(see)
library(performance)
### check assumptions
check_model(m2)
check_model(m4)

library(multcompView)
library(emmeans)
library(multcomp)
# estimated marginal means
mm4 = emmeans(m4, ~ treatment*variety, adjust = "none", type = "response")
mm4

plot(mm4, comparisons=T)
pvalues <- pairs(mm4, adjust="none")

cld(mm4,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", ### no adjustment = LSD for multiple comparisons
    method = "paiwise") 


#calculate summary statistics of variables 
library(plyr)
flwr.stats1 <- ddply(flowering, c("treatment", 'variety', "dat", 'vartreat'),
                     summarise,
                     N    = length(Plants.Flowering),#We use length instead of count.
                     sum = sum(Plants.Flowering))
flwr.stats1

levels(flwr.stats1$variety)[levels(flwr.stats1$variety)=='Tahoe'] <- 'Tahoe Cinco'

ggplot(flwr.stats1, aes(x=dat, y=sum, fill = treatment)) + 
  facet_wrap(~ variety) +
  geom_bar(stat = 'identity', position="dodge")+ 
  scale_fill_manual(values = c("#1B9E77", '#E7298A'))+
  #geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, size=1,position=position_dodge(0)) +
  theme_classic() +
  ylab(bquote('Plants flowering')) +
  xlab('Treatment') +
  labs(color = 'Days after transplanting (DAT)', linetype = 'Cultivar') + #rename legend title
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title=element_text(size=14)) + 
  theme(strip.text = element_text(size = 14)) + 
  theme(axis.title.y = element_text(vjust=0.9)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 10), breaks = 0:10*2)+
  scale_x_discrete(expand = expansion(add=c(0.2,0.2)))+
  #geom_text(aes(label=letters.biomass, y = mean + se), vjust = -.5, size = 6)+
  theme(legend.position = 'right')

bar.flwr <- ggplot(flwr.stats1, aes(x=dat, y=sum, fill = treatment)) + 
  facet_wrap(~ variety) +
  geom_bar(stat = 'identity', position="dodge", color = "black")+ 
  scale_fill_grey()+  
  theme_classic() +
  ylab('Plants flowering') +
  xlab('Days after transplanting (DAT)') +
  labs(fill = 'Treatment') + #rename legend title
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 18)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  theme(legend.text=element_text(size=18)) + 
  theme(legend.title=element_text(size=18)) + 
  theme(strip.text = element_text(size = 18)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 9), breaks = 0:10*2)
  #theme(legend.position = 'none') +
  #geom_text(aes(label=l, y = sum + 1), vjust = -.5, size = 7) 
bar.flwr




########### Cannabinoids - 10% and 90% pistil die back ##############
#### CBD #####
read.csv('10-90 pistil die back.csv')
pistildieback <- read.csv('10-90 pistil die back.csv', header = T, sep = ",")

#make dates into a date format  
pistildieback$date <- as.Date(pistildieback$Harvest.Date,'%m/%d/%Y')

#create DAT column 
pistildieback$dat <- as.factor(difftime(as.POSIXct(pistildieback$date), as.POSIXct('2021-06-14', tz="PT"), units="days"))

###Check to see that all combinations are correct
xtabs(~ variety + dat, data = pistildieback)
xtabs(~ treatment + dat, data = pistildieback)
xtabs(~ plot + dat, data = pistildieback)
xtabs(~ treatment + variety, data = pistildieback)
xtabs(~ block + variety, data = pistildieback)

# as.factor
pistildieback$variety <- as.factor(paste(pistildieback$variety))
pistildieback$treatment <- as.factor(paste(pistildieback$treatment))

#Combined Variety & Treatment 
pistildieback$vartreat <- as.factor(paste(pistildieback$variety, pistildieback$treatment))

# combined dat - die.back
pistildieback$percdat <- as.factor(paste(pistildieback$die.back, pistildieback$dat))

pistildieback$CBD <- pistildieback$X.CBDA.as.CBD
pistildieback$THC <- pistildieback$X.THCA.as.THC

# #identify outliers 
# library(mvoutlier)
# x <- weekly.fresh$vartreat
# y <- weekly.fresh$X.CBDA
# boxplot(y~x, xlab='Variety', ylab='Canopy Cover (m)')
# stripchart(y~x, 
#            vertical=TRUE,
#            method='jitter',
#            pch=21, col='blue',
#            bg='red',
#            add=TRUE)
# #outliers <- identify(y~x)

# outliers
pistildieback4 <- pistildieback[-c(25, 36, 1, 12),]
pistildieback4

library(viridis)
# boxplot
ggplot(pistildieback, aes(x=dat, y=CBD, color = vartreat)) + 
  geom_violin() +
  scale_fill_viridis(discrete = F, alpha=0.6) +
  geom_jitter(size=3, alpha=0.9) +
  ylab('% CBD') +
  xlab('Days After Transplanting') +
  labs(color = 'Treatment', shape = 'Treatment') + #rename legend title
  theme(axis.title.x = element_text(size = 14), axis.text = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  ggtitle('Cannabinoids - Field 2021') +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=14)) +
  theme(strip.text = element_text(size = 14))

ggplot(pistildieback, aes(x=dat, y=THC, color = vartreat)) + 
  geom_violin() +
  scale_fill_viridis(discrete = F, alpha=0.6) +
  geom_jitter(size=3, alpha=0.9) +
  ylab('% THC') +
  xlab('Days After Transplanting') +
  labs(color = 'Treatment', shape = 'Treatment') + #rename legend title
  theme(axis.title.x = element_text(size = 14), axis.text = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  ggtitle('Cannabinoids - Field 2021') +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=14)) +
  theme(strip.text = element_text(size = 14))

#calculate summary statistics of variables 
library(plyr)
cbd.perc <- ddply(pistildieback4, c("vartreat", "dat", "die.back"),
                  summarise,
                  N    = length(CBD),#We use length instead of count.
                  mean = mean(CBD),
                  sd   = sd(CBD),
                  se   = sd / sqrt(N))
cbd.perc

cbd.perc.var.dieback <- ddply(pistildieback4, c("vartreat", "die.back", "treatment", "variety"),
                              summarise,
                              N    = length(CBD),#We use length instead of count.
                              mean = mean(CBD),
                              sd   = sd(CBD),
                              se   = sd / sqrt(N))
cbd.perc.var.dieback
#write.table(cbd.perc, file = 'CBD pistil die back field 2021 Means.csv', sep = ",", quote = FALSE, row.names = F)

# create labels 
cbd.perc.var.dieback <- cbd.perc.var.dieback %>% mutate(l = case_when(
  die.back == "90%" & vartreat == "Red Bordeaux Control" ~ "c", 
  die.back == "90%" & vartreat == "Red Bordeaux N+" ~ "bc",
  die.back == "90%" & vartreat == "Berry Blossom N+" ~ "bc",
  die.back == "90%" & vartreat == "Berry Blossom Control" ~ "bc",
  die.back == "90%" & vartreat == "Tahoe Cinco Control" ~ "ab",
  die.back == "90%" & vartreat == "Tahoe Cinco N+" ~ "a"))

line.perc.cbd.bw <- ggplot(cbd.perc.var.dieback, aes(x=die.back, y=mean, group = vartreat, shape = vartreat, linetype = vartreat)) + 
  geom_point(size = 4) +
  geom_line(size = 1) +
  scale_shape_manual(labels = c('BBc', 'BBn', 'RBc', 'RBn', 'TCc', 'TCn'), values=c(15, 0, 17, 2, 18, 9))+
  scale_linetype_discrete(labels = c('BBc', 'BBn', 'RBc', 'RBn', 'TCc', 'TCn')) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.08, size=1) +
  theme_classic() +
  ylab('CBD (% dry weight)') +
  xlab('Pistil dieback') +
  theme(axis.title.y = element_text(vjust=0.9)) +
  theme(axis.title.x = element_text(size = 16), axis.text = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=16)) + 
  theme(strip.text = element_text(size = 16)) +
  #scale_y_continuous(expand = c(0,0))+
  #theme(legend.position = "none") +
  #geom_text_repel(aes(label=l), size = 7, nudge_y  = subset(cbd.perc.var.dieback, die.back == '10%')$mean) +
  geom_text_repel(aes(label=l), size = 7, nudge_x = 0.2) +
  theme(legend.position = c(.1, .8)) + 
  theme(legend.box.background = element_rect(color="black", size=1))
line.perc.cbd.bw

ggplot(cbd.perc.var.dieback, aes(x=die.back, y=mean, group = vartreat, color = treatment, linetype = variety)) + 
  geom_point(size = 4) +
  geom_line(size = 1) +
  scale_linetype_manual(values=c('solid', "dotted", "dashed"))+
  scale_color_manual(values = c("#1B9E77", '#E7298A'))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, size=1,position=position_dodge(0)) +
  theme_classic() +
  ylab('CBD (% dry weight)') +
  xlab('Pistil dieback') +
  labs(color = 'Treatment', linetype = 'Cultivar') + #rename legend title
  theme(axis.title.x = element_text(size = 16), axis.text = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  theme(legend.text=element_text(size=16)) + 
  #theme(legend.title = element_blank()) +
  theme(legend.title=element_text(size=14)) +
  theme(strip.text = element_text(size = 16)) +
  scale_y_continuous(expand = c(0,0),limits = c(7, 14))+
  geom_text_repel(aes(label=l), size = 7, nudge_x = 0.2)+
  #geom_text_repel(aes(x=die.back, y=mean, label=l), size = 7, nudge_x = c(-0.2, 0.2))+
  theme(legend.position = c(.1, .8)) + 
  theme(legend.box.background = element_rect(color="black", size=1))+
  theme(legend.position = 'none')

# remove days 10%
cbd.perc.var.dieback.90 <- subset(cbd.perc.var.dieback, die.back != '10%')
cbd.perc.var.dieback.90

letters.cbd = c("bc", "bc", "c", "bc", "ab", "a")

bar.cbd.90 <- ggplot(cbd.perc.var.dieback.90, aes(x = treatment, y=mean, fill = treatment)) +
  facet_wrap(~ variety) +
  #facet_grid(vars(Leaf.Position), vars(variety), scales = "free", space = 'free') + 
  geom_bar(stat = 'identity', position="dodge", color = "black")+ 
  scale_fill_grey()+
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab("CBD (% dry weight)") + 
  xlab("Treatment") +
  ggtitle('90% pistil dieback harvest') + 
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 18)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  theme(legend.text=element_text(size=18)) + 
  theme(legend.title=element_text(size=18)) + 
  theme(strip.text = element_text(size = 18)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,15)) + 
  coord_cartesian(ylim = c(8, 15))+
  theme(legend.position = 'none') +
  geom_text(aes(label=letters.cbd, y = mean + se), vjust = -.5, size = 7)  
bar.cbd.90

# manuscript CBD figure 
#ggsave(bar.cbd.90, file="field22 cbd.90.bw.TIFF",width=13, height=6,dpi=600,path="H:/WNM Conference 2023")

### CBD pistil die back models ###
library(lme4)
p0 <- lmer(CBD ~ die.back*treatment + die.back*variety + die.back*variety + (1|plot), data= pistildieback4, REML = F )
p1 <- lmer(CBD ~ die.back*treatment + die.back*variety + (1|plot), data = pistildieback4, REML = F)
#p2 <- lmer(CBD ~ die.back*treatment + die.back*variety + (1|block), data = pistildieback4, REML = F)
#p3 <- lmer(CBD ~ die.back*treatment*variety + (1|plot), data = pistildieback4, REML = F)
#p5 <- lmer(CBD ~ die.back*treatment*variety + (1|block) , data = pistildieback4, REML = F)
#p6 <- lmer(CBD ~ die.back*treatment + die.back*variety +  (1|plot)+ (1|dat), data = pistildieback4, REML = F)
#p7 <- lmer(CBD ~ die.back*treatment + die.back*variety + (1|block) + (1|dat), data = pistildieback4, REML = F)
#p3.1 <- lmer(CBD ~ die.back*treatment + die.back*variety + (1|plot) , data = pistildieback4, REML = F)
#p6.1 <- lmer(CBD ~ die.back*treatment + die.back*variety +  (1|block), data = pistildieback4, REML = F)

check_model(p0)
AIC(p0, p1)
check_model(p1)
#check_model(p3)
anova(p0, p1)

#p8 <- lmer(sqrt(CBD) ~ die.back*treatment*variety + (1|block), data = pistildieback4, REML = F)
#p9 <- lmer(sqrt(CBD) ~ die.back*treatment*variety + (1|plot) , data = pistildieback4, REML = F)
#p10 <- lmer(sqrt(CBD) ~ die.back*treatment + die.back*variety + (1|plot) + (1|dat), data = pistildieback4, REML = F)
#p11 <- lmer(sqrt(CBD) ~ die.back*treatment + die.back*variety + (1|plot) + (1|dat), data = pistildieback4, REML = F)
#p12 <- lmer(sqrt(CBD) ~ die.back*treatment + die.back*variety + (1|block)  + (1|dat), data = pistildieback4, REML = F)
p10.1 <- lmer(sqrt(CBD) ~ die.back*treatment + die.back*variety + treatment*variety + (1|plot) , data = pistildieback4, REML = F)
#p11.1 <- lmer(sqrt(CBD) ~ die.back*treatment + die.back*variety + (1|block) , data = pistildieback4, REML = F)
AIC(p9, p10.1)
check_model(p10.1)

#p14 <- lmer(log(CBD) ~ die.back*treatment*variety + (1|block) , data = pistildieback4, REML = F)
#p15 <- lmer(log(CBD) ~ die.back*treatment + die.back*variety + (1|block) + (1|dat), data = pistildieback4, REML = F)
p13.1 <- lmer(log(CBD) ~ die.back*treatment + die.back*variety + treatment*variety + (1|plot), data = pistildieback4, REML = F)
#p14.1 <- lmer(log(CBD) ~ die.back*treatment + die.back*variety + (1|block), data = pistildieback4, REML = F)
check_model(p13.1)

library(rcompanion)
cbd.transfrom = transformTukey(pistildieback4$CBD, plotit=T)

cbd.sqr = (-1)*(pistildieback4$CBD)^(0.3)   # Avoid complex numbers
p0 <- lmer(cbd.sqr ~ die.back*treatment + treatment*variety + die.back*variety + (1|plot), data= pistildieback4, REML = F )
#p1 <- lmer(cbd.sqr ~ die.back*treatment + treatment*variety + die.back*variety + (1|block), data= pistildieback4, REML = F )
p2 <- lmer(cbd.sqr ~ die.back*treatment + die.back*variety + (1|plot), data = pistildieback4, REML = F)
#p3 <- lmer(cbd.sqr ~ die.back*treatment + die.back*variety + (1|block), data = pistildieback4, REML = F)
p4 <- lmer(cbd.sqr ~ die.back*treatment*variety + (1|plot), data = pistildieback4, REML = F)
#p5 <- lmer(cbd.sqr ~ die.back*treatment*variety + (1|block) , data = pistildieback4, REML = F)
#p6 <- lmer(cbd.sqr ~ die.back*treatment + die.back*variety +  (1|plot)+ (1|dat), data = pistildieback4, REML = F)
#p7 <- lmer(cbd.sqr ~ die.back*treatment + die.back*variety + (1|block) + (1|dat), data = pistildieback4, REML = F)

AIC(p0, p2, p4)

library(merTools)
RMSE.merMod(p0)
RMSE.merMod(p2)

library(see)
library(performance)
### check assumptions
check_model(p0)
check_model(p2)

Anova(p0)
Anova(p2)

RMSE.merMod(p10.1)
RMSE.merMod(p13.1)

library(see)
library(performance)
### check assumptions
check_model(p13.1)
check_model(p10.1)
check_model(p1)

library(phia)
plot(interactionMeans(p10.1))
plot(interactionMeans(p2))
plot(interactionMeans(p13.1))

#Anova(p13.1)
Anova(p10.1, by = "die.back")
#Anova(p1)
#Anova(p0)

library(multcompView)
library(emmeans)
library(multcomp)
# estimated marginal means
e4 = emmeans(p2, specs = pairwise ~ variety*treatment|die.back,  adjust = "none")
e4
e3 = emmeans(p10.1, ~ die.back*treatment + die.back*variety + treatment*variety, adjust = "none")
e3
e1 = emmeans(p10.1, specs = pairwise ~ variety*treatment|die.back, adjust = "none", type = "response") # compare levels of treatment within group
e1

cld(e1,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")    ### no adjustment = LSD for multiple comparisons

# cld comparisons for both 10% and 90%
cld(e3,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    #type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")    ### no adjustment = LSD for multiple comparisons

contrast(e1, method="pairwise",adjust="none", type = "response", infer = T)
pairs(e3, alpha = 0.05, adjust = "none") # contrast stats without CI's

pwpp(e1, by = "die.back", adjust = "none")
pairs(e3, adjust = "none")
plot(e1, comparisons = T, adjust = "none", by = "die.back") # if red arrows overlap groups are not sig from each other




####### THC ########## 
thc.perc <- ddply(pistildieback4, c("vartreat", "dat", "die.back", "treatment", "variety"),
                  summarise,
                  N    = length(THC),#We use length instead of count.
                  mean = mean(THC),
                  sd   = sd(THC),
                  se   = sd / sqrt(N))
thc.perc
#write.table(thc.perc, file = 'THC pistil die back field 2021 Means.csv', sep = ",", quote = FALSE, row.names = F)

ggplot(thc.perc, aes(x=die.back, y=mean, group = vartreat, color = treatment, linetype = variety)) + 
  geom_point(size = 4) +
  geom_line(size = 1) +
  scale_linetype_manual(values=c('solid', "dotted", "dashed"))+
  scale_color_manual(values = c("#1B9E77", '#E7298A'))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, size=1,position=position_dodge(0)) +
  theme_classic() +
  ylab('THC (% dry weight)') +
  xlab('Pistil dieback') +
  labs(color = 'Treatment', linetype = 'Cultivar') + #rename legend title
  theme(axis.title.x = element_text(size = 16), axis.text = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  theme(legend.text=element_text(size=16)) + 
  #theme(legend.title = element_blank()) +
  geom_hline(yintercept = 0.3, linetype = "dashed", size = 1)+
  theme(legend.title=element_text(size=14)) +
  theme(strip.text = element_text(size = 16)) +
  scale_y_continuous(expand = c(0,0),limits = c(0.25, 0.50))+
  #geom_text_repel(aes(label=l), size = 7, nudge_x = 0.2)+
  #geom_text_repel(aes(x=die.back, y=mean, label=l), size = 7, nudge_x = c(-0.2, 0.2))+
  theme(legend.position = c(.1, .8)) + 
  theme(legend.box.background = element_rect(color="black", size=1))+
  theme(legend.position = 'none')

bar.thc <- ggplot(thc.perc, aes(factor(x = die.back), y=mean, fill = vartreat)) +
  #facet_grid(~ die.back, scales = 'free_x' , space = "free_x") +
  geom_col(position = position_dodge2(preserve = 'single')) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge2(preserve = 'single')) +
  #scale_fill_manual(values=c('blue','darkred', 'green', 'red', 'darkorange', 'black', '#56B4E9')) +
  theme_classic() +
  ylab('THC %' ) +
  xlab('Days after transplanting') +
  labs(fill = 'Variety') + #rename legend title
  theme(axis.title.x = element_text(size = 14), axis.text = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14))  +
  #ggtitle('Cannabinoids Week 7 & 9 - Greenhouse Hemp 2021') +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=14)) +
  theme(axis.title.x = element_text(size = 14), axis.text = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.55)) +
  #geom_text(aes(label=thc.l, y = mean + se), vjust = -.5, size = 6) + 
  theme(legend.position = "none") 
#scale_x_discrete(expand = expansion(add=c(0.2,0.2)))
bar.thc

library(lme4)
library(rcompanion)
thc.transfrom = transformTukey(pistildieback4$THC, plotit=T)

thc.sqr = (pistildieback4$CBD)^(0.875)   # Avoid complex numbers

t0 <- lmer(thc.sqr ~ die.back*treatment + treatment*variety + die.back*variety + (1|plot), data= pistildieback4, REML = F )
#p1 <- lmer(thc.sqr ~ die.back*treatment + treatment*variety + die.back*variety + (1|block), data= pistildieback4, REML = F )
t2 <- lmer(thc.sqr ~ die.back*treatment + die.back*variety + (1|plot), data = pistildieback4, REML = F)
#p3 <- lmer(thc.sqr ~ die.back*treatment + die.back*variety + (1|block), data = pistildieback4, REML = F)
#t4 <- lmer(thc.sqr ~ die.back*treatment*variety + (1|plot), data = pistildieback4, REML = F)
#p5 <- lmer(thc.sqr ~ die.back*treatment*variety + (1|block) , data = pistildieback4, REML = F)
#p6 <- lmer(thc.sqr ~ die.back*treatment + die.back*variety +  (1|plot)+ (1|dat), data = pistildieback4, REML = F)
#p7 <- lmer(thc.sqr ~ die.back*treatment + die.back*variety + (1|block) + (1|dat), data = pistildieback4, REML = F)
AIC(t0, t2, t4)
Anova(t0)
Anova(t2)
library(see)
library(performance)
### check assumptions
check_model(t0)
check_model(t2)


pt0 <- lmer(THC ~ die.back*treatment + treatment*variety + die.back*variety + (1|plot), data = pistildieback4, REML = F)
pt1 <- lmer(THC ~ die.back*treatment + die.back*variety + (1|plot), data = pistildieback4, REML = F)
pt2 <- lmer(THC ~ die.back*treatment + treatment*variety + die.back*variety + (1|block), data = pistildieback4, REML = F)
pt3 <- lmer(THC ~ die.back*treatment + die.back*variety + (1|block) , data = pistildieback4, REML = F)
AIC(pt0, pt1, pt2, pt3)

check_model(pt0)
check_model(pt1) #-120
check_model(pt3) # -120

pt7 <- lmer(sqrt(THC) ~ die.back*treatment + die.back*variety + (1|plot) , data = pistildieback4, REML = F)
pt10 <- lmer(sqrt(THC) ~ die.back*treatment+die.back*variety + treatment*variety + (1|plot) , data = pistildieback4, REML = F)
AIC(pt7,pt10)
check_model(pt7)
check_model(pt10)

#pt13 <- lmer(log(THC) ~ die.back*treatment + die.back*variety + (1|block) , data = pistildieback4, REML = F)
pt14 <- lmer(log(THC) ~ die.back*treatment + die.back*variety + (1|plot) , data = pistildieback4, REML = F)
pt15 <- lmer(log(THC) ~ die.back*treatment + die.back*variety + treatment*variety+ (1|plot) , data = pistildieback4, REML = F)
#pt16 <- lmer(log(THC) ~ die.back*treatment*variety + (1|block) , data = pistildieback4, REML = F)
#pt17 <- lmer(log(THC) ~ die.back*treatment*variety + (1|plot) , data = pistildieback4, REML = F)
AIC(pt14, pt15)

check_model(pt14)
check_model(pt15)

plot(interactionMeans(pt2))

library(merTools)
RMSE.merMod(pt0)
RMSE.merMod(pt1)
RMSE.merMod(pt2)
RMSE.merMod(pt3)
RMSE.merMod(pt7) #low
RMSE.merMod(pt10) #low
RMSE.merMod(pt14)
RMSE.merMod(pt15)


anova(pt7, pt10)

#Anova(pt10)
#Anova(pt15)
Anova(pt7)


# estimated marginal means comparisons 
e1 = emmeans(t2, ~ die.back*treatment + die.back*variety, adjust = "none", type = "response") # compare levels of treatment within group
e1
e2 = emmeans(pt10, specs = pairwise ~ variety*treatment|die.back, adjust = "none", type = "response")
e2
e3 = emmeans(pt2, ~ die.back+treatment + die.back*variety, adjust = "none")
e3
et4 = emmeans(t0, specs = pairwise ~ variety*treatment|die.back,  adjust = "none")
et4

# 10% and 90% die back comparisons 
cld(e1,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")    ### no adjustment = LSD for multiple comparisons

# seperate 10% and 90% pistil die back
cld(et4,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")    ### no adjustment = LSD for multiple comparisons


# confidence intervals with statistical tests 
contrast(e2, method="pairwise",adjust="none", type = "response", infer = T)
pairs(e3, alpha = 0.05, adjust = "none") # contrast stats without CI's
#write.csv(SPAD.contrasts, file = 'SPAD.contrasts.csv', quote = FALSE, row.names = T)

pwpp(e3, by = "treatment", adjust = "none")
pairs(e3, adjust = "none")
plot(e2, comparisons = T, adjust = "none", by = "variety") # if red arrows overlap groups are not sig from each other

# create labels 
thc.perc <- thc.perc %>% mutate(l = case_when(
  die.back == "90%" & vartreat == "Red Bordeaux Control" ~ "d", 
  die.back == "90%" & vartreat == "Red Bordeaux N+" ~ "cd",
  die.back == "90%" & vartreat == "Berry Blossom Control" ~ "cd",
  die.back == "90%" & vartreat == "Berry Blossom N+" ~ "bc",
  die.back == "90%" & vartreat == "Tahoe Cinco Control" ~ "ab",
  die.back == "90%" & vartreat == "Tahoe Cinco N+" ~ "a"))

line.perc.thc.bw <- ggplot(thc.perc, aes(x=die.back, y=mean, group = vartreat, shape = vartreat, linetype = vartreat)) + 
  geom_point(size = 4) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.08, size=1) +
  theme_classic() +
  ylab('THC (% dry weight)') +
  xlab('Pistil dieback') +
  geom_hline(yintercept=0.3, linetype='dashed', col = 'black')+
  scale_shape_manual(labels = c('BBc', 'BBn', 'RBc', 'RBn', 'TCc', 'TCn'), values=c(15, 0, 17, 2, 18, 9))+
  scale_linetype_discrete(labels = c('BBc', 'BBn', 'RBc', 'RBn', 'TCc', 'TCn')) +
  theme(axis.title.x = element_text(size = 16), axis.text = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  theme(legend.title = element_blank()) +
  theme(legend.text=element_text(size=16)) + 
  #theme(legend.title=element_text(size=14)) +
  theme(strip.text = element_text(size = 16)) +
  #scale_y_continuous(expand = c(0,0))+
  #theme(legend.position = "none") +
  geom_text_repel(aes(label=l), size = 7, nudge_x = 0.2) +
  theme(legend.position = c(.1, .8)) + 
  theme(legend.box.background = element_rect(color="black", size=1))
line.perc.thc.bw



######### merge CBD & THC 10% and 90% die back manuscript figures #############
library(gridExtra)
grid.arrange(line.perc.cbd, line.perc.thc, ncol = 2, nrow = 1)
cbd.thc.bw <- grid.arrange(line.perc.cbd.bw, line.perc.thc.bw, ncol = 2, nrow = 1)
cbd.thc.bw

# figure with labels 
cbd.thc.bw3 <- plot_grid(line.perc.cbd.bw, line.perc.thc.bw, labels = "AUTO", label_size = 18)
cbd.thc.bw3

#ggsave(cbd.thc.bw3, file="cbd.thc.bw3.TIFF",width=12, height=6,dpi=600,path="C:/Users/mfarnisa/OneDrive - University of Nevada, Reno/Thesis/Figures")


############### CBD-to-THC pistil die back Ratios #########
ratio <- ddply(pistildieback4, c("vartreat", "die.back", "treatment", "variety"),
               summarise,
               N    = length(Ratio),#We use length instead of count.
               mean = mean(Ratio),
               sd   = sd(Ratio),
               se   = sd / sqrt(N))
ratio

ggplot(ratio, aes(x=die.back, y=mean, group = vartreat, color = treatment, linetype = variety)) + 
  geom_point(size = 4) +
  geom_line(size = 1) +
  scale_linetype_manual(values=c('solid', "dotted", "dashed"))+
  scale_color_manual(values = c("#1B9E77", '#E7298A'))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, size=1,position=position_dodge(0)) +
  theme_classic() +
  ylab('CBD:THC') +
  xlab('Pistil dieback') +
  labs(color = 'Treatment', linetype = 'Cultivar') + #rename legend title
  theme(axis.title.x = element_text(size = 16), axis.text = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  theme(legend.text=element_text(size=16)) + 
  #theme(legend.title = element_blank()) +
  theme(legend.title=element_text(size=14)) +
  theme(strip.text = element_text(size = 16)) +
  scale_y_continuous(expand = c(0,0),limits = c(22, 32))+
  #geom_text_repel(aes(label=l), size = 7, nudge_x = 0.2)+
  #geom_text_repel(aes(x=die.back, y=mean, label=l), size = 7, nudge_x = c(-0.2, 0.2))+
  theme(legend.position = c(.1, .8)) + 
  theme(legend.box.background = element_rect(color="black", size=1))+
  theme(legend.position = 'none')

library(lme4)
library(rcompanion)

r.transfrom = transformTukey(pistildieback4$Ratio, plotit=T)

r.sqr = (-1)*(pistildieback4$CBD)^(-2.7)   # Avoid complex numbers

r0 <- lmer(r.sqr ~ die.back*treatment + treatment*variety + die.back*variety + (1|plot), data= pistildieback4, REML = F )
#r1 <- lmer(r.sqr ~ die.back*treatment + treatment*variety + die.back*variety + (1|block), data= pistildieback4, REML = F )
r2 <- lmer(r.sqr ~ die.back*treatment + die.back*variety + (1|plot), data = pistildieback4, REML = F)
#r3 <- lmer(r.sqr ~ die.back*treatment + die.back*variety + (1|block), data = pistildieback4, REML = F)
r4 <- lmer(r.sqr ~ die.back*treatment*variety + (1|plot), data = pistildieback4, REML = F)
#r5 <- lmer(r.sqr ~ die.back*treatment*variety + (1|block) , data = pistildieback4, REML = F)
#r6 <- lmer(r.sqr ~ die.back*treatment + die.back*variety +  (1|plot)+ (1|dat), data = pistildieback4, REML = F)
#r7 <- lmer(r.sqr ~ die.back*treatment + die.back*variety + (1|block) + (1|dat), data = pistildieback4, REML = F)
AIC(r0, r2, r4)
Anova(r0)
Anova(r2)
library(see)
library(performance)
### check assumptions
check_model(r0)
check_model(r2)  

r1 <- lmer(Ratio ~ die.back*treatment + treatment*variety + die.back*variety + (1|block), data = pistildieback4, REML = F)
AIC(r1)

r2 <- lmer(log(Ratio) ~ die.back*treatment + die.back*variety + (1|block), data = pistildieback4, REML = F)
r3 <- lmer(sqrt(Ratio) ~ die.back*treatment + die.back*variety + (1|block), data = pistildieback4, REML = F)

check_model(r1)

check_model(r2)
check_model(r3)

Anova(r1)
plot(interactionMeans(r1))


rr1 = emmeans(r2, ~ die.back*treatment + treatment*variety + die.back*variety, adjust = "none")
rr1
rr2 = emmeans(r2, specs = pairwise ~ variety*treatment|die.back, adjust = "none", type = "response") # compare levels of treatment within group
rr2

cld(rr2,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")    ### no adjustment = LSD for multiple comparisons

# create labels 
ratio <- ratio %>% mutate(l = case_when(
  die.back == "10%" & vartreat == "Red Bordeaux Control" ~ "x", 
  die.back == "10%" & vartreat == "Berry Blossom N+" ~ "z",
  die.back == "10%" & vartreat == "Red Bordeaux N+" ~ "xy",
  die.back == "10%" & vartreat == "Tahoe Cinco Control" ~ "yz",
  die.back == "10%" & vartreat == "Berry Blossom Control" ~ "xy",
  die.back == "10%" & vartreat == "Tahoe Cinco N+" ~ "z",
  die.back == "90%" & vartreat == "Red Bordeaux Control" ~ "a", 
  die.back == "90%" & vartreat == "Red Bordeaux N+" ~ "a",
  die.back == "90%" & vartreat == "Berry Blossom Control" ~ "a",
  die.back == "90%" & vartreat == "Berry Blossom N+" ~ "b",
  die.back == "90%" & vartreat == "Tahoe Cinco Control" ~ "b",
  die.back == "90%" & vartreat == "Tahoe Cinco N+" ~ "b"))

#write.table(ratio, file = 'Ratio.CBD.THC pistil die back field 2021 Means.csv', sep = ",", quote = FALSE, row.names = F)
line.ratio.bw <- ggplot(ratio, aes(x=die.back, y=mean, group = vartreat, shape = vartreat, linetype = vartreat)) + 
  geom_point(size = 4) +
  geom_line(size = 1) +
  scale_shape_manual(labels = c('BBc', 'BBn', 'RBc', 'RBn', 'TCc', 'TCn'), values=c(15, 0, 17, 2, 18, 9))+
  scale_linetype_discrete(labels = c('BBc', 'BBn', 'RBc', 'RBn', 'TCc', 'TCn')) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.4, size=1,position=position_dodge(0)) +
  theme_classic() +
  scale_fill_grey()+
  ylab('CBD-to-THC ratio') +
  #ylab(bquote('Plant height  '(cm~plant^-1))) +
  xlab('Pistil dieback') +
  #labs(shape = 'Cultivar - Treatment', linetype = 'Cultivar - Treatment') + #rename legend title
  theme(axis.title.x = element_text(size = 16), axis.text = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  #ggtitle('Plant Height linear model - Field 2021') +
  theme(legend.text=element_text(size=16)) + 
  theme(legend.title = element_blank()) +
  #theme(legend.title=element_text(size=14)) +
  theme(strip.text = element_text(size = 16)) +
  scale_y_continuous(expand = c(0,0),limits = c(22.5, 34))+
  geom_text_repel(aes(x=die.back, y=mean, label=l), size = 7, nudge_x = c(-0.2, 0.2))+
  theme(legend.position = c(.1, .8)) + 
  theme(legend.box.background = element_rect(color="black", size=1))
#theme(legend.position = 'none')
line.ratio.bw

# remove days 10%
ratio.90 <- subset(ratio, die.back != '10%')
ratio.90

letters.ratio = c("a", "b", "a", "a", "b", "b")

bar.ratio.90 <- ggplot(ratio.90, aes(x = treatment, y=mean, fill = treatment)) +
  facet_wrap(~ variety) +
  #facet_grid(vars(Leaf.Position), vars(variety), scales = "free", space = 'free') + 
  geom_bar(stat = 'identity', position="dodge", color = "black")+ 
  scale_fill_grey()+
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab("CBD-to-THC Ratio") + 
  xlab("Treatment") +
  ggtitle('90% pistil dieback harvest') + 
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 18)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  theme(legend.text=element_text(size=18)) + 
  theme(legend.title=element_text(size=18)) + 
  theme(strip.text = element_text(size = 18)) + 
  scale_y_continuous(expand = c(0,0),limits = c(0, 35)) +
  coord_cartesian(ylim = c(25, 33))+
  theme(legend.position = 'none') +
  geom_text(aes(label=letters.ratio, y = mean + se), vjust = -.5, size = 7)  
bar.ratio.90

#ggsave(bar.ratio.90, file="field22 ratio.90.bw.TIFF",width=13, height=6,dpi=600,path="H:/WNM Conference 2023")




############# CBD yield/biomass ##################
read.csv('CBD.Biomass.Yield.csv')
cbd.biomass <- read.csv('CBD.Biomass.Yield.csv', header = T, sep = ",")

#make dates into a date format  
cbd.biomass$date <- as.Date(cbd.biomass$Harvest.Date,'%m/%d/%Y')

#create DAT column 
cbd.biomass$dat <- as.factor(difftime(as.POSIXct(cbd.biomass$date), as.POSIXct('2021-06-14', tz="PT"), units="days"))

#factor variety
cbd.biomass$variety <- as.factor(cbd.biomass$variety)

#factor treatment 
cbd.biomass$treatment <- as.factor(cbd.biomass$treatment)

# factor row
cbd.biomass$row <- as.factor(cbd.biomass$row)
cbd.biomass$block <- as.factor(cbd.biomass$block)

cbd.biomass$vartreat <- as.factor(paste(cbd.biomass$variety, cbd.biomass$treatment))

# remove outliers 1 & 12 like pistildieback4
cbd.biomass4 <- cbd.biomass[-c(1, 12),]

xtabs(~ treatment + variety, data = cbd.biomass4)

# calculate biomass in g/plant
cbd.biomass4$g.plant = (cbd.biomass4$Flower.Total.Dried_g/2)*(cbd.biomass4$X.CBDA.as.CBD/100)
#cbd.biomass4$g = (cbd.biomass4$Flower.Total.Dried_g/2)*(cbd.biomass4$X.CBDA.as.CBD)
cbd.biomass4$kg = ((cbd.biomass4$Flower.Total.Dried_g/2)*(cbd.biomass4$X.CBDA.as.CBD))/1000

#calculate summary statistics of variables 
flower.cbd <- ddply(cbd.biomass4, c("variety", "treatment",'vartreat'),
                    summarise,
                    N    = length(g.plant),#We use length instead of count.
                    mean = mean(g.plant),
                    sd   = sd(g.plant),
                    se   = sd / sqrt(N))
flower.cbd
#write.table(flower.cbd, file = 'CBD.biomass yield 90% pistil die back field 2021 Means.csv', sep = ",", quote = FALSE, row.names = F)

ggplot(flower.cbd, aes(x=treatment, y=mean, fill = treatment)) + 
  facet_wrap(~ variety) +
  geom_bar(stat = 'identity', position="dodge")+ 
  scale_fill_manual(values = c("#1B9E77", '#E7298A'))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, size=1,position=position_dodge(0)) +
  theme_classic() +
  ylab(bquote('CBD yield  '(g~plant^-1))) + 
  xlab('Treatment') +
  labs(color = 'Treatment', linetype = 'Cultivar') + #rename legend title
  theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(legend.text=element_text(size=12)) + 
  theme(legend.title=element_text(size=14)) + 
  theme(strip.text = element_text(size = 14)) + 
  theme(axis.title.y = element_text(vjust=0.9)) +
  scale_y_continuous(breaks=0:75*10, expand = c(0,0), limits = c(0,65)) + 
  scale_x_discrete(expand = expansion(add=c(0.2,0.2)))+
  #geom_text(aes(label=cbd.biomass.l, y = mean + se), vjust = -.5, size = 7)+
  theme(legend.position = 'none')


interaction.plot(x.factor     = cbd.biomass4$treatment,
                 trace.factor = cbd.biomass4$variety,
                 response     = cbd.biomass4$g.plant,
                 fun = mean,
                 type="b",
                 col=c("black","red","green"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

library(lme4)
g1 <- lmer(g.plant ~ treatment*variety + (1|block) , data = cbd.biomass4, REML = F)
#g2 <- lmer(g.plant ~ treatment+variety + (1|block) , data = cbd.biomass4, REML = F)
g3 <- lmer(sqrt(g.plant) ~ treatment*variety + (1|block), data = cbd.biomass4, REML = F)
#g4 <- lmer(sqrt(g.plant) ~ treatment+variety + (1|block) , data = cbd.biomass4, REML = F)
g5 <- lmer(log(g.plant) ~ treatment*variety + (1|block) , data = cbd.biomass4, REML = F)
#g6 <- lmer(log(g.plant) ~ treatment+variety + (1|block) , data = cbd.biomass4, REML = F)

AIC(g1, g2)
# g1  8 168.1760
# g2  6 176.5996

plot(interactionMeans(g1))

library(see)
library(performance)
### check assumptions
check_model(g1)# ok?
#check_model(g3)# no
check_model(g5) #no 

library(merTools)
RMSE.merMod(g1)
RMSE.merMod(g3)
RMSE.merMod(g5)

#Anova(g1)
Anova(g3)
#Anova(g5)

anova(g1, g5)

# estimated marginal means
e2 = emmeans(g5, ~ treatment*variety, adjust = "none")
e2
e3 = emmeans(g3, ~ treatment*variety, adjust = "none")
e3

cld(e2,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")    ### no adjustment = LSD for multiple comparisons

cld(e3,
    alpha   = 0.05,
    reversed  = T, ### reverse the order of letters
    Letters = letters,    ### Use lower-case letters for .group
    type    = "response", ### Report emmeans in orginal scale
    adjust =  "none", 
    method = "pairwise")    ### no adjustment = LSD for multiple comparisons

cbd.biomass.l = c("c", "b", "c", "bc", "c", "a")

bar.cbd.yield.bw <- ggplot(flower.cbd, aes(x = treatment, y=mean, fill = treatment))+
  facet_wrap(~ variety) +
  geom_bar(stat = 'identity', position="dodge", color = "black")+ 
  scale_fill_grey()+
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, position=position_dodge(0.9))+
  theme_classic() + 
  ylab(bquote('CBD yield  '(g~plant^-1))) +
  #xlab("N Treatment") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(vjust=0.9)) +
  theme(legend.text=element_text(size=16)) + 
  theme(legend.title=element_text(size=16)) + 
  theme(strip.text = element_text(size = 16)) + 
  #scale_y_continuous(breaks=0:70*10, expand = c(0,0), limits = c(0,70)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 65)) +
  theme(legend.position = 'none') +
  geom_text(aes(label=cbd.biomass.l, y = mean + se), vjust = -.5, size = 7) 
bar.cbd.yield.bw


##### merge CBD yield and CBD-to-THC Ratio manuscript figures ####
yield.ratio.bw2 <- plot_grid(bar.cbd.yield.bw, line.ratio.bw, labels = "AUTO", label_size = 18)
yield.ratio.bw2
#ggsave(yield.ratio.bw2, file="yield.ratio.bw.TIFF",width=12, height=6,dpi=600,path="C:/Users/mfarnisa/OneDrive - University of Nevada, Reno/Thesis/Figures")
