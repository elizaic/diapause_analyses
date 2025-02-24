##############################################################################.
############.
############         Estimate location adaption of D. carinulata
############      diapasuse timing in two environmnets
############.      
############                Written by: ##############
############            Date Last Modified: 2/24/2025
############.      
##############################################################################.



# Packages ----------------------------------------------------------------

library(readxl)
library(tidyverse)
library(lubridate)
library(survival)
library(survminer)
library(emmeans)
library(ggfortify)
library(viridis)
library(glmmTMB)
library(DHARMa)
library(car)
library(grid)
library(ggeffects)
library(broom)
# library(RColorBrewer)  

select <- dplyr::select


# Data ------------------------------------------------------------------------


## Diapause timing -----------------------------------------------------------
# load diapause data, delete missing trials, 
# create time column, make trt & pop factors

diapauseData <- read_csv("data/diapause_local_adapt_data.csv")
diapauseData <- diapauseData %>% mutate(
  start_date = mdy(start_date),
  end_date = mdy(end_date)
)
str(diapauseData)
diapauseData2 <- diapauseData %>%
  filter(diapause != "missing") %>%
  filter(is.na(notes) == T) %>%
  select(-notes)  %>% 
  mutate(time = (end_date - start_date)) %>%
  mutate(treatment = as_factor(treatment)) %>%
  mutate(population = as_factor(population)) %>%
  mutate(diapause = as.numeric(diapause)) %>%
  mutate(core_edge = case_when(
    population == "Bl" | population == "Lc" | population == "Lj" | population == "Wi" ~ "edge",
    population == "De" | population == "Hu" | population == "Lo" | population == "Pu" ~ "core"
  )) %>%
  mutate(core_edge = as_factor(core_edge)) %>%
  mutate(time = as.numeric(time, units = 'days'))

  
## Weight of beetles ---------------------------------------------------------

#load weight data
weightData <- read_csv("data/weight_local_adapt_data.csv")
weightData <- weightData %>% mutate(
  eclosion = mdy(eclosion))

#Combine data sets
diapauseData3 <- left_join(diapauseData2, 
          weightData %>% select(beetle_ID, weight), 
                 by = "beetle_ID") %>%
  mutate(latitude = case_when(
    population == "Bl" ~ 33.912,
    population == "Lc" ~ 34.593,
    population == "Lj" ~ 34.342,
    population == "Wi" ~ 34.422,
    population == "De" ~ 39.144,
    population == "Hu" ~ 40.063,
    population == "Lo" ~ 44.856,
    population == "Pu" ~ 38.268
  )) %>%
  mutate(elevation = case_when(
    population == "Bl" ~ 96.9,
    population == "Lc" ~ 1669.56,
    population == "Lj" ~ 1433,
    population == "Wi" ~ 1204.74,
    population == "De" ~ 1386.43,
    population == "Hu" ~ 1189.62,
    population == "Lo" ~ 1115.35,
    population == "Pu" ~ 1448.19
  )) %>%
  mutate(degree_days = case_when(
    population == "Bl" ~ 4524,
    population == "Lc" ~ 2043,
    population == "Lj" ~ 2547,
    population == "Wi" ~ 3516,
    population == "De" ~ 1767,
    population == "Hu" ~ 1841,
    population == "Lo" ~ 1438,
    population == "Pu" ~ 1990
  )) %>%
  mutate(distFromCore = 39.144 - latitude, 
         distFromEdge = 34.422 - latitude,
         time = ifelse(diapause == 0 & time >= 41, 43, time))
diapauseData3$population <- factor(diapauseData3$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))
diapauseData3$core_edge <- factor(diapauseData3$core_edge, levels = c('core', 'edge'))


str(diapauseData3)

#Find range of ages when beetles started treatment
ageAtStart <- left_join(diapauseData2 %>% select(beetle_ID, start_date),
          weightData %>% select(beetle_ID, eclosion), 
          by = "beetle_ID") 
range(ageAtStart$start_date - ageAtStart$eclosion, na.rm = T)

## weight of LJ population (checking for high weight beetles that might indicated non-carinulata)
diapauseData3 %>% filter(population == 'Lj') %>% 
  ggplot(aes(x = weight)) +
  geom_histogram(bins = 50)


## Climate and site characteristics --------------------------------------------
site_char <- data.frame(
  population = c("Bl", "Lc", "Lj", "Wi", "De", "Hu", "Lo", "Pu"),
  latitude = c(33.912,34.593,34.342,34.422,39.144,40.063,44.856,38.268),
  elevation = c(96.9, 1669.56, 1433, 1204.74,1386.43, 1189.62, 1115.35, 1448.19),
  degree_days = c(4524, 2043, 2547, 3516, 1767, 1841, 1438, 1990),
  core_edge = c(rep("1", 4), rep("0", 4)),
  letter = c("H", "E", "G", "F", "C", "B", "A", "D")
  )

site_char$core_edge <- as.numeric(site_char$core_edge)
str(site_char)


# Objects for plots, etc.

edge_col = viridis(3)[1]
core_col = viridis(3)[2]

# Analyses ------------------------------------------------------------------

## Days to diapause ---------------------------------------------------------

### Diapausers only ---------------------------------------------------------

## dataset
diapauseDataDiapausing <- diapauseData3 %>%
  filter(diapause == 1)

## summary statistics
dtd_n2 <- diapauseDataDiapausing %>% group_by(treatment, population) %>% summarise(
  mean = mean(time),
  sd = sd(time),
  n = n()
)
dtd_n2 <- dtd_n2 %>% mutate(n2 = replace(n, n==20, 'n=20') )

diapauseDataDiapausing %>% group_by(treatment, core_edge) %>% summarise(
  mean = mean(time),
  sd = sd(time),
  n = n()
)
dtd_n1 <- c("n = 90", '15', "79", '126')

diapauseData3%>% group_by(treatment, core_edge) %>% summarise(
  mean = mean(time),
  sd = sd(time),
  n = n()
)


## model with core_Edge, pop as rand
modDiapausing <- glmmTMB(time ~ core_edge * treatment + (1| population), data = diapauseDataDiapausing,
                         family = poisson)
simulateResiduals(modDiapausing) %>% plot()
summary(modDiapausing)
Anova(modDiapausing, type = 3)
emmeans(modDiapausing, pairwise ~ core_edge|treatment, type = 'response')
emmeans(modDiapausing, pairwise ~ treatment|core_edge, type = 'response')

# estimated marginal means for core and edge
emoutDiapausing <- emmeans(modDiapausing, pairwise ~ core_edge*treatment, type = "response")
meansDiapausing_df <- as.data.frame(emoutDiapausing$emmeans)
meansDiapausing_df$core_edge <- factor(meansDiapausing_df$core_edge, levels = c('core', 'edge'))
meansDiapausing_df <- cbind(meansDiapausing_df,dtd_n1)

diapauseDataDiapausing$core_edge <- factor(diapauseDataDiapausing$core_edge, levels = c('core', 'edge'))


## model with populations - low sample sizes here!
modDiapausingpop <- glmmTMB(time ~  treatment * population, data = diapauseDataDiapausing,
                            family = poisson)
simulateResiduals(modDiapausingpop) %>% plot()
summary(modDiapausingpop)
Anova(modDiapausingpop, type = 3)
emmeans(modDiapausingpop, pairwise ~ treatment*population, type = 'response')

# estimated marginal means for each population
emoutDiapausingPops <- emmeans(modDiapausingpop, pairwise ~ population*treatment, type = "response")
meansDiapausingPops_df <- as.data.frame(emoutDiapausingPops$emmeans) %>%
  mutate(core_edge = case_when(
    population == "Bl" | population == "Lc" | population == "Lj" | population == "Wi" ~ "edge",
    population == "De" | population == "Hu" | population == "Lo" | population == "Pu" ~ "core"
  ))
meansDiapausingPops_df$population <- factor(meansDiapausingPops_df$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))
diapauseDataDiapausing$population <- factor(diapauseDataDiapausing$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))



## Plot of core/edge and population means
# Figure S1
ggplot(data = meansDiapausing_df, aes(x = treatment, y = rate, fill = core_edge)) +
  #add data points
  geom_point(data = diapauseDataDiapausing, aes(x = treatment, y = time, color = core_edge),
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2), alpha = 0.5) +
  
  #population lines and points
  geom_line(data = meansDiapausingPops_df, aes(group = population),
            position = position_dodge(width = 0.5), alpha = 0.5) +
  geom_pointrange(data = meansDiapausingPops_df, aes(x = treatment, y = rate, ymin = asymp.LCL, ymax = asymp.UCL, group = population, fill = core_edge),
                  position = position_dodge(width = 0.5),
                  color = 'black', alpha = .5, shape = 22, size = 0.75) +
  
  #core/edge lines and points
  geom_line(aes(group = core_edge), position = position_dodge(width = 0.5), linewidth = 1.2) +
  geom_pointrange(aes(x = treatment, y = rate, ymin = asymp.LCL, ymax = asymp.UCL),
                  position = position_dodge(width = 0.5), shape = 21, size = 1) +
  
  #sample size
  geom_text(aes(label = dtd_n1, y = -1.5), position = position_dodge(width = 0.5)) +
  
  #formatting
  labs(x = "Diapause-Inducing Daylength", y = "Days to Diapause (95% CI)", fill = "Collection \nLocation") +
  scale_x_discrete(labels = c('long' = "North", 'short' = "South")) +
  scale_color_manual(values = c('core' = 'grey75', 'edge' = 'grey25'), labels = c('edge' = "South", 'core' = "North")) +
  guides(color = 'none') +
  scale_fill_manual(values = c('core' = 'white', 'edge' = 'black'), labels = c('edge' = "South", 'core' = "North")) +
  theme_bw(base_size = 15)
ggsave("plots/S1_days_only_diapausers.png", dpi = 400, width = 7, height = 8)





### All samples --------------------------------------------------------------

## summary statistics
summ_dtd_all_pop <- diapauseData3 %>% group_by(treatment, population) %>% summarise(
  mean = mean(time),
  sd = sd(time),
  n = n()
)

mean(summ_dtd_all_pop$n)

summ_dtd_all_pop <- summ_dtd_all_pop %>% mutate(n2 = replace(n, n==26, 'n=26') )

summ_dtd_all <- diapauseData3 %>% group_by(treatment, core_edge) %>% summarise(
  mean = mean(time),
  sd = sd(time),
  n = n()
)
summ_dtd_all <- summ_dtd_all %>% mutate(n2 = replace(n, n==170, 'n = 170         ') )



## model with core_Edge, pop as rand
modAllSamp <- glmmTMB(time ~ core_edge * treatment + (1| population), data = diapauseData3,
                      family = poisson())
simulateResiduals(modAllSamp) %>% plot()
summary(modAllSamp)
Anova(modAllSamp, type = 3)
emmeans(modAllSamp, pairwise ~ core_edge|treatment, type = 'response')
emmeans(modAllSamp, pairwise ~ treatment|core_edge, type = 'response')


emoutmodAllSamp <- emmeans(modAllSamp, pairwise ~ core_edge*treatment, type = "response")
emoutmodAllSamp_df <- as.data.frame(emoutmodAllSamp$emmeans)
emoutmodAllSamp_df$core_edge <- factor(emoutmodAllSamp_df$core_edge, levels = c('core', 'edge'))



# Population analyses - all samples
modAllSamppop <- glmmTMB(time ~  treatment * population, data = diapauseData3,
                         family = poisson())
simulateResiduals(modAllSamppop) %>% plot()
summary(modAllSamppop)
Anova(modAllSamppop, type = 3)
emmeans(modAllSamppop, pairwise ~ treatment|population, type = 'response')
emmeans(modAllSamppop, pairwise ~ population|treatment, type = 'response')

emoutAllSampPops <- emmeans(modAllSamppop, pairwise ~ population*treatment, type = "response")
emoutAllSampPops_df <- as.data.frame(emoutAllSampPops$emmeans) %>%
  mutate(core_edge = case_when(
    population == "Bl" | population == "Lc" | population == "Lj" | population == "Wi" ~ "edge",
    population == "De" | population == "Hu" | population == "Lo" | population == "Pu" ~ "core"
  ),
  letter = case_when(
    population == "Bl" ~ "H",
    population == "Lc" ~ "E",
    population == "Lj" ~ "G",
    population == "Wi" ~ "F",
    population == "De" ~ "C",
    population == "Hu" ~ "B",
    population == "Lo" ~ "A",
    population == "Pu" ~ "D"
  ))
emoutAllSampPops_df$population <- factor(emoutAllSampPops_df$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))


## Combining population and core edge plots in reaction norm style
days_plot1 <- ggplot(data = emoutmodAllSamp_df, aes(x = treatment, y = rate, fill = core_edge)) +
  #box for non-diapausers
  annotate('rect', xmin = c(0.7,1.7), xmax = c(1.3, 2.3), ymin = 42, ymax = 44, fill = 'grey95', color = "white") +
  # 
  # #add data points
  geom_point(data = diapauseData3, aes(x = treatment, y = time, color = core_edge),
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2), alpha = 0.5) +
  
  #population lines and points
  geom_line(data = emoutAllSampPops_df, aes(group = population),
            position = position_dodge(width = 0.5), alpha = 0.5) +
  geom_pointrange(data = emoutAllSampPops_df, aes(x = treatment, y = rate, ymin = asymp.LCL, ymax = asymp.UCL, group = population),
                  position = position_dodge(width = 0.5),
                  color = 'grey50', fill = 'white', shape = 22, size = 1.5) +
  
  #add letters for populations
  geom_text(data = emoutAllSampPops_df, aes(label = letter, y = rate, group = population), 
            position = position_dodge(width = 0.5)) +
  
  #core/edge lines and points
  geom_line(aes(group = core_edge), position = position_dodge(width = 0.5), linewidth = 1.2) +
  geom_pointrange(aes(x = treatment, y = rate, ymin = asymp.LCL, ymax = asymp.UCL),
                  position = position_dodge(width = 0.5), shape = 21, size = 1.5) +
  
  #sample size
  geom_text(data = summ_dtd_all, aes(label = n2, y = -1.5), position = position_dodge(width = 0.35)) +
  
  #non-diapauser label
  annotate("label", x = c(1,2), y = 45.5, label = "Non-Diapausers", fill = 'grey95', label.size = NA) +
  
  #formatting
  labs(x = "Diapause-Inducing Daylength", y = "Days to Diapause (95% CI)", fill = "Collection \nLocation") +
  scale_x_discrete(labels = c('long' = "North", 'short' = "South")) +
  scale_color_manual(values = c('core' = 'grey75', 'edge' = 'grey25'), labels = c('edge' = "South", 'core' = "North")) +
  guides(color = 'none') +
  scale_fill_manual(values = c('core' = 'white', 'edge' = 'black'), labels = c('edge' = "Edge/South", 'core' = "Core/North")) +
  ylim(-1.5,52)+
  
  theme_bw(base_size = 15)
days_plot1





### Survival analysis -------------------------------------------------------


## Create survival object
survObj <- Surv(diapauseData3$time, diapauseData3$diapause)

## Population
model1c <- survreg(survObj ~ population * treatment, dist = "gaussian", data = diapauseData3)
model1c
anova(model1c)
emmeans(model1c, pairwise ~ treatment | population, type = "response")
survreg_means <- emmeans(model1c, pairwise ~ treatment | population, type = "response")
survreg_means_df <- as.data.frame(survreg_means$emmeans) %>%
  mutate(core_edge = case_when(
    population == "Bl" | population == "Lc" | population == "Lj" | population == "Wi" ~ "edge",
    population == "De" | population == "Hu" | population == "Lo" | population == "Pu" ~ "core"
  )) %>%
  unite("treatment_core_edge", c(treatment,core_edge), remove = F)
survreg_means_df$population <- factor(survreg_means_df$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))


## Core and Edge
model1b_2 <- survreg(survObj ~ core_edge * treatment + population, dist ="gaussian", data = diapauseData3)
model1b_2
summary(model1b_2)
emmeans(model1b_2, pairwise ~ treatment | core_edge, type = "response")
survreg_means_CE <- emmeans(model1b_2, pairwise ~ treatment | core_edge, type = "response")
survreg_means_CE_df <- as.data.frame(survreg_means_CE$emmeans) %>%
  unite("treatment_core_edge", c(treatment,core_edge), remove = F)
survreg_means_CE_df$core_edge <- factor(survreg_means_CE_df$core_edge, levels = c('core', 'edge'))



## Combining population and core edge plots in reaction norm style
ggplot(data = survreg_means_CE_df, aes(x = treatment, y = emmean, fill = core_edge)) +
  #box for non-diapausers
  annotate('rect', xmin = c(0.7,1.7), xmax = c(1.3, 2.3), ymin = 42, ymax = 44, fill = 'grey95', color = "white") +
  
  #add data points
  geom_point(data = diapauseData3, aes(x = treatment, y = time, color = core_edge),
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2), alpha = 0.5) +
  
  #population lines and points
  geom_line(data = survreg_means_df, aes(group = population),
            position = position_dodge(width = 0.5), alpha = 0.5) +
  geom_pointrange(data = survreg_means_df, aes(x = treatment, y = emmean, ymin = lower.CL, ymax = upper.CL, group = population, fill = core_edge),
                  position = position_dodge(width = 0.5),
                  color = 'black', alpha = .5, shape = 22, size = 0.75) +
  
  #core/edge lines and points
  geom_line(aes(group = core_edge), position = position_dodge(width = 0.5), size = 1.2) +
  geom_pointrange(aes(x = treatment, y = emmean, ymin = lower.CL, ymax = upper.CL),
                  position = position_dodge(width = 0.5), shape = 21, size = 1) +
  
  #non-diapauser label
  annotate("label", x = c(1,2), y = 47, label = "Non-Diapausers", fill = 'grey95', label.size = NA) +
  
  #sample size
  # geom_text(aes(label = dtd_n1, y = -1.5), position = position_dodge(width = 0.5)) +
  
  #formatting
  labs(x = "Diapause-Inducing Daylength", y = "Days to Diapause (95% CI)", fill = "Collection \nLocation") +
  scale_x_discrete(labels = c('long' = "North", 'short' = "South")) +
  scale_color_manual(values = c('core' = 'grey75', 'edge' = 'grey25'), labels = c('edge' = "South", 'core' = "North")) +
  guides(color = 'none') +
  scale_fill_manual(values = c('core' = 'white', 'edge' = 'black'), labels = c('edge' = "South", 'core' = "North")) +
  theme_bw(base_size = 15)

ggsave("plots/S2_days_survival.png", dpi = 400, width = 7, height = 8)





## Proportion in diapause ----------------------------------------------------


# Summary statistics
diapusePropSummaryPop <- diapauseData3 %>% group_by(population, treatment) %>%
  summarise(
    n = n(),
    prop = mean(diapause),
    sd = sd(diapause),
    se = sd/sqrt(n),
    latitude = mean(latitude)
  )
diapusePropSummaryPop

diapauseData3 %>% group_by(treatment, core_edge) %>% summarise(
  mean_time = mean(time),
  prop = mean(diapause),
  sd_prop = sd(diapause),
  n = n()
)

##Sample size 
prop_sample_size <- diapauseData3 %>% group_by(population, treatment) %>% summarize (
  n = n()
)
mean(prop_sample_size$n)

diapauseData3$population <- factor(diapauseData3$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))
diapauseData3$core_edge <- factor(diapauseData3$core_edge, levels = c('core', 'edge'))



# Core/edge Logistic regression

#Model core/edge
logmodel1 <- glmmTMB(diapause ~ core_edge * treatment + (1|population), family = "binomial", data = diapauseData3)
simulateResiduals(logmodel1) %>% plot()
summary(logmodel1)
Anova(logmodel1, type = 3)
# ranef(logmodel1, condVar = F)
# coef(logmodel1)
emmeans(logmodel1, pairwise ~ core_edge | treatment, type = "response")


logistic_means_ce <- emmeans(logmodel1, pairwise ~ treatment * core_edge, type = "response")
logistic_means_ce_df <- as.data.frame(logistic_means_ce$emmeans) %>%
  unite("treatment_core_edge", c(treatment,core_edge), remove = F)

 

#Population models

## Bayesian approach to fix complete separation
logmodel_bayes <- bayesglm(diapause ~ population * treatment, family = 'binomial', data = diapauseData3)
vif(logmodel_bayes)

summary(logmodel_bayes)
Anova(logmodel_bayes, type = 3)
emmeans(logmodel_bayes, pairwise ~ treatment | population, type = "response")
emmeans(logmodel_bayes, pairwise ~ population | treatment, type = "response")

logisticbayes_means_pop <- emmeans(logmodel_bayes, pairwise ~ treatment * population, type = "response")
logisticbayes_means_pop_df <- as.data.frame(logisticbayes_means_pop$emmeans) %>%
  mutate(core_edge = case_when(
    population == "Bl" | population == "Lc" | population == "Lj" | population == "Wi" ~ "edge",
    population == "De" | population == "Hu" | population == "Lo" | population == "Pu" ~ "core"
  ))
logisticbayes_means_pop_df$population <- factor(logisticbayes_means_pop_df$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))

logisticbayes_means_pop_df <- logisticbayes_means_pop_df %>% mutate(letter = rep(c("A", "B", "C", "D", "E", "F", "G", "H"), each = 2))


## Combining population and core edge plots in reaction norm style
prop_plot1 <- ggplot(data = logistic_means_ce_df, aes(x = treatment, y = prob, fill = core_edge)) +
  
  #population lines and points
  geom_line(data = logisticbayes_means_pop_df, aes(x = treatment, y = prob, group = population),
            position = position_dodge(width = 0.5), alpha = 0.5) +
  geom_pointrange(data = logisticbayes_means_pop_df, aes(x = treatment, y = prob, ymin = asymp.LCL, ymax = asymp.UCL, group = population, fill = core_edge),
                  position = position_dodge(width = 0.5),
                  color = 'grey50', fill = "white", shape = 22, size = 1.5) +
  
  #add letters for populations
  geom_text(data = logisticbayes_means_pop_df, aes(label = letter, y = prob, group = population),
            position = position_dodge(width = 0.5)) +
  
  #core/edge lines and points
  geom_line(aes(group = core_edge), position = position_dodge(width = 0.5), size = 1.2) +
  geom_pointrange(aes(x = treatment, y = prob, ymin = asymp.LCL, ymax = asymp.UCL),
                  position = position_dodge(width = 0.5), shape = 21, size = 1.5) +
  
  #sample size
  # geom_text(aes(label = dtd_n1, y = -1.5), position = position_dodge(width = 0.5)) +
  
  #formatting
  labs(x = "Diapause-Inducing Daylength", y = "Proportion in Diapause (95% CI)", fill = "Collection \nLocation") +
  scale_x_discrete(labels = c('long' = "North", 'short' = "South")) +
  scale_color_manual(values = c('core' = 'grey75', 'edge' = 'grey25'), labels = c('edge' = "South", 'core' = "North")) +
  guides(color = 'none') +
  scale_fill_manual(values = c('core' = 'white', 'edge' = 'black'), labels = c('edge' = "Edge/South", 'core' = "Core/North")) +
  theme_bw(base_size = 15)
prop_plot1

#Combine proportion in diapause with days until diapause plots
ggarrange(prop_plot1, days_plot1, nrow = 1, labels = "AUTO", common.legend = T, legend = 'bottom')

ggsave("plots/Fig2_prop_days.png", dpi = 400, width = 12, height = 8)






##Correlation between prop & days ---------------------------------------------


##Using population means

logisticbayes_means_pop_df

emoutAllSampPops_df

prop_days_pop <- full_join(logisticbayes_means_pop_df,emoutAllSampPops_df, by = c('treatment', 'population'))
cor.test(prop_days_pop$prob, prop_days_pop$rate, method = c("pearson"))

ggplot(data = prop_days_pop)+
  geom_smooth(aes(x = prob, y = rate), method = lm, color = 'black') +
  geom_jitter(aes(x = prob, y = rate, fill = core_edge.x, shape = treatment),
              size = 3, width = 0.02) +
  # geom_text(aes(label = letter, x = prob+.02, y = rate + 1), size = 4) +
  labs(y = "Days until Diapause", x = "Proportion in Diapause", fill = 'Origin',
       shape = 'Treatment') +
  annotate("text", y = 40, x = 0.75, label = "italic(R) ^ 2 == -0.990", parse = TRUE, size = 6) +
  annotate("text", y = 37, x = 0.75, label = "italic(P) < 0.001", parse = TRUE, size = 5) +
  # annotate("text", y = -2, x = site_char$latitude[1:4], label = site_char$letter[1:4], size = 4) +
  scale_shape_manual(values = c(21, 22), labels = c('long' = "North/Long day", 'short' = "South/Short day")) +
  scale_fill_manual(values = c('core' = 'white', 'edge' = 'black'), labels = c('edge' = "Edge/South", 'core' = "Core/North")) +
  guides(fill = guide_legend('Origin', override.aes = list(shape = 21))) +
  theme_bw(base_size = 15)

ggsave("plots/S3_corr_days_prop.png", dpi = 400, width = 7, height = 8)


prop_days_pop_long <- prop_days_pop %>% filter(treatment == 'long')
cor.test(prop_days_pop_long$prob , prop_days_pop_long$rate , method = c("pearson"))

prop_days_pop_short <- prop_days_pop %>% filter(treatment == 'short')
cor.test(prop_days_pop_short$prob , prop_days_pop_short$rate , method = c("pearson"))



## Environmental Variable regressions ----------------------------------------

##Are latitude, elevation, and degree days correlated at the 8 sites?

#EDGE
cor.test(site_char$latitude[5:8], site_char$elevation[5:8])
cor.test(site_char$latitude[5:8], site_char$degree_days[5:8])
cor.test(site_char$elevation[5:8], site_char$degree_days[5:8])

#CORE
cor.test(site_char$latitude[1:4], site_char$elevation[1:4])
cor.test(site_char$latitude[1:4], site_char$degree_days[1:4])
cor.test(site_char$elevation[1:4], site_char$degree_days[1:4])

#For all sites: Latitude and elevation are not correlated among the 8 sites (Pearson's=0.131, P=0.758). Degree days is negatively correlated (Pearson's) with both latitude (-0.716, P=0.046) and elevation (-0.743, P=0.035).

#For edge sites: Latitude and elevation are not sig. correlated (-0.861, P=0.14). Degree days is correlated with latitude (-0.953, P=0.047), but not altitude (0.771, P=0.229)

#For core sites: Latitude and elevation are correlated (0.964, P=0.04). Degree days is not sig. correlated with latitude (-0.879, P=0.12) nor altitude (-0.948, P= 0.052)

e_s <- diapauseData3 %>% filter(core_edge == 'edge') %>% filter(treatment == 'short')


mod_e_s_lat <- lm(time ~ latitude, data = e_s)
summary(mod_e_s_lat)

mod_e_s_ele <- lm(time ~ elevation, data = e_s)
summary(mod_e_s_ele)

mod_e_s_deg <- lm(time ~ degree_days, data = e_s)
summary(mod_e_s_deg)
AIC(mod_e_s_lat,mod_e_s_ele,mod_e_s_deg)



##Plots edge in short treatment


e_lat <- ggplot(data = e_s, aes(x = latitude, y = time)) +
  geom_smooth(method = lm, color = 'black') +
  geom_jitter(width = .01, height = 0.2) +
  labs(y = "Days to Diapause", x = "Latitude") +
  annotate("text", y = 37, x = 34.2, label = "italic(R) ^ 2 == 0.226", parse = TRUE, size = 6) +
  annotate("text", y = -2, x = site_char$latitude[1:4], label = site_char$letter[1:4], size = 4) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())
e_lat

e_ele <- ggplot(data = e_s, aes(x = elevation, y = time)) +
  geom_smooth(method = lm, color = 'black') +
  geom_jitter(width = 15, height = 0.2) +
  theme_bw(base_size = 15) +
  labs(y = "Days to Diapause", x = "Elevation") +
  annotate("text", y = 37, x = 850, label = "italic(R) ^ 2 == 0.340", parse = TRUE, size = 6) +
  annotate("text", y = -2, x = site_char$elevation[1:4], label = site_char$letter[1:4], size = 4) +
  theme(panel.grid = element_blank())
e_ele

e_deg <- ggplot(data = e_s, aes(x = degree_days, y = time)) +
  geom_smooth(method = lm, color = 'black') +
  geom_jitter(width = 15, height = 0.2) +
  theme_bw(base_size = 15) +
  labs(y = "Days to Diapause", x = "Cumulative Degree Days") +
  annotate("text", y = 37, x = 2750, label = "italic(R) ^ 2 == 0.494", parse = TRUE, size = 6) +
  annotate("text", y = -2, x = site_char$degree_days[1:4], label = site_char$letter[1:4], size = 4) +
  theme(panel.grid = element_blank())
e_deg

e_plots <- ggarrange(e_lat,e_ele,e_deg, nrow = 1) %>% annotate_figure(top = text_grob("Edge Sites in Southern Environment", face = "bold", size = 15))
#export: 1750 x 700 px



## Core, Long 
c_l <- diapauseData3 %>% filter(core_edge == 'core') %>% filter(treatment == 'long')

mod_c_l_all1 <- lm(time ~ latitude+degree_days+elevation, data = c_l)
vif(mod_c_l_all1)

mod_c_l_all <- lm(time ~ latitude+degree_days, data = c_l)
summary(mod_c_l_all)
Anova(mod_c_l_all, type = 3)
plot(simulateResiduals(mod_c_l_all))




mod_c_l_lat <- lm(time ~ latitude, data = c_l)
summary(mod_c_l_lat)

mod_c_l_ele <- lm(time ~ elevation, data = c_l)
summary(mod_c_l_ele)

mod_c_l_deg <- lm(time ~ degree_days, data = c_l)
summary(mod_c_l_deg)

AIC(mod_c_l_lat,mod_c_l_ele,mod_c_l_deg)


##Plots core in long treatment

c_lat <- ggplot(data = c_l, aes(x = latitude, y = time)) +
  geom_smooth(method = lm, color = 'black') +
  geom_jitter(width = .1, height = 0.2) +
  theme_bw(base_size = 15) +
  labs(y = "Days to Diapause", x = "Latitude") +
  annotate("text", y = 37, x = 41, label = "italic(R) ^ 2 == 0.071", parse = TRUE, size = 6) +
  annotate("text", y = -2, x = site_char$latitude[5:8], label = site_char$letter[5:8], size = 4) +
  theme(panel.grid = element_blank())
c_lat

c_ele <- ggplot(data = c_l, aes(x = elevation, y = time)) +
  geom_smooth(method = lm, color = 'black') +
  geom_jitter(width = 5, height = 0.2) +
  theme_bw(base_size = 15) +
  labs(y = "Days to Diapause", x = "Elevation") +
  annotate("text", y = 37, x = 1250, label = "italic(R) ^ 2 == 0.084", parse = TRUE, size = 6) +
  annotate("text", y = -2, x = site_char$elevation[5:8], label = site_char$letter[5:8], size = 4) +
  theme(panel.grid = element_blank())
c_ele

c_deg <- ggplot(data = c_l, aes(x = degree_days, y = time)) +
  geom_smooth(method = lm, color = 'black') +
  geom_jitter(width = 10, height = 0.2) +
  theme_bw(base_size = 15) +
  labs(y = "Days to Diapause", x = "Cumulative Degree Days") +
  annotate("text", y = 37, x = 1600, label = "italic(R) ^ 2 == 0.036", parse = TRUE, size = 6) +
  annotate("text", y = -2, x = site_char$degree_days[5:8], label = site_char$letter[5:8], size = 4) +
  theme(panel.grid = element_blank())
c_deg

c_plots <- ggarrange(c_lat,c_ele,c_deg, nrow = 1) %>% annotate_figure(top = text_grob("Core Sites in Northern Environment", face = "bold", size = 15))

ggarrange(c_plots, e_plots, nrow = 2, labels = 'AUTO')

ggsave("plots/Fig3_enviro_in_home.png", dpi = 400, width = 12, height = 10)




## Long all 
long_all <- diapauseData3 %>% filter(treatment == 'long')
#check vif
mod_long_all1 <- lm(time ~ latitude+degree_days+elevation, data = long_all)
vif(mod_long_all1)

mod_long_all_lat <- lm(time ~ latitude, data = long_all)
summary(mod_long_all_lat)

mod_long_all_ele <- lm(time ~ elevation, data = long_all)
summary(mod_long_all_ele)

mod_long_all_deg <- lm(time ~ degree_days, data = long_all)
summary(mod_long_all_deg)



## Short all 
short_all <- diapauseData3 %>% filter(treatment == 'short')
#check vif
mod_short_all1 <- lm(time ~ latitude+degree_days+elevation, data = short_all)
vif(mod_short_all1)

mod_short_all_lat <- lm(time ~ latitude, data = short_all)
summary(mod_short_all_lat)

mod_short_all_ele <- lm(time ~ elevation, data = short_all)
summary(mod_short_all_ele)

mod_short_all_deg <- lm(time ~ degree_days, data = short_all)
summary(mod_short_all_deg)


##Plots in long treatment

long_lat <- ggplot(data = long_all, aes(x = latitude, y = time)) +
  geom_smooth(method = lm, color = 'black') +
  geom_jitter(width = .1, height = 0.2) +
  theme_bw(base_size = 15) +
  labs(y = "Days to Diapause", x = "Latitude") +
  annotate("text", y = 37, x = 41, label = "italic(R) ^ 2 == 0.179", parse = TRUE, size = 6) +
  annotate("text", y = -2, x = site_char$latitude, label = site_char$letter, size = 4) +
  theme(panel.grid = element_blank())
long_lat

long_ele <- ggplot(data = long_all, aes(x = elevation, y = time)) +
  geom_smooth(method = lm, color = 'black') +
  geom_jitter(width = 5, height = 0.2) +
  theme_bw(base_size = 15) +
  labs(y = "Days to Diapause", x = "Elevation") +
  annotate("text", y = 37, x = 750, label = "italic(R) ^ 2 == -0.004", parse = TRUE, size = 6) +
  annotate("text", y = -2, x = site_char$elevation, label = site_char$letter, size = 4) +
  theme(panel.grid = element_blank())
long_ele

long_deg <- ggplot(data = long_all, aes(x = degree_days, y = time)) +
  geom_smooth(method = lm, color = 'black') +
  geom_jitter(width = 10, height = 0.2) +
  theme_bw(base_size = 15) +
  labs(y = "Days to Diapause", x = "Cumulative Degree Days") +
  annotate("text", y = 25, x = 3250, label = "italic(R) ^ 2 == 0.092", parse = TRUE, size = 6) +
  annotate("text", y = -2, x = site_char$degree_days, label = site_char$letter, size = 4) +
  theme(panel.grid = element_blank())
long_deg


long_plots <- ggarrange(long_lat,long_ele,long_deg, nrow = 1) %>% annotate_figure(top = text_grob("All Sites in Northern Environment", face = "bold", size = 15))
#export: 1750 x 700 px



##Plots edge in short treatment


short_lat <- ggplot(data = short_all, aes(x = latitude, y = time)) +
  geom_smooth(method = lm, color = 'black') +
  geom_jitter(width = .01, height = 0.2) +
  theme_bw(base_size = 15) +
  labs(y = "Days to Diapause", x = "Latitude") +
  annotate("text", y = 37, x = 38, label = "italic(R) ^ 2 == 0.139", parse = TRUE, size = 6) +
  annotate("text", y = -2, x = site_char$latitude, label = site_char$letter, size = 4) +
  theme(panel.grid = element_blank())
short_lat

short_ele <- ggplot(data = short_all, aes(x = elevation, y = time)) +
  geom_smooth(method = lm, color = 'black') +
  geom_jitter(width = 15, height = 0.2) +
  theme_bw(base_size = 15) +
  labs(y = "Days to Diapause", x = "Elevation") +
  annotate("text", y = 37, x = 850, label = "italic(R) ^ 2 == 0.326", parse = TRUE, size = 6) +
  annotate("text", y = -2, x = site_char$elevation, label = site_char$letter, size = 4) +
  theme(panel.grid = element_blank())
short_ele

short_deg <- ggplot(data = short_all, aes(x = degree_days, y = time)) +
  geom_smooth(method = lm, color = 'black') +
  geom_jitter(width = 15, height = 0.2) +
  theme_bw(base_size = 15) +
  labs(y = "Days to Diapause", x = "Cumulative Degree Days") +
  annotate("text", y = 37, x = 2750, label = "italic(R) ^ 2 == 0.546", parse = TRUE, size = 6) +
  annotate("text", y = -2, x = site_char$degree_days, label = site_char$letter, size = 4) +
  theme(panel.grid = element_blank())
short_deg



short_plots <- ggarrange(short_lat,short_ele,short_deg, nrow = 1) %>% annotate_figure(top = text_grob("All Sites in Southern Environment", face = "bold", size = 15))

ggarrange(long_plots, short_plots, nrow = 2, labels = 'AUTO')

ggsave("plots/S4_enviro_all.png", dpi = 400, width = 12, height = 10)



