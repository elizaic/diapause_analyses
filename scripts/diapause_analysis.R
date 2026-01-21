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
library(arm)
# library(RColorBrewer) 
library(ggpattern)
library(tinytable)


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
  first_frost = c(346, 286, 302, 322, 281, 273, 275, 288),
  core_edge = c(rep("1", 4), rep("0", 4)),
  letter = c("H", "E", "G", "F", "C", "B", "A", "D")
  )

site_char$core_edge <- as.numeric(site_char$core_edge)
str(site_char)


# Analyses ------------------------------------------------------------------

# Function to extract p-values from Anova table and format for plotting
plot_pval_labs <- function(model, 
                           term_labels = c("Core/Edge", "Treatment", "Core/Edge*Treatment")) {
  labs <- Anova(model, type = 3) %>%
    as.data.frame() %>%
    mutate( 
      rowname = rownames(.),
      pval = `Pr(>Chisq)`,
      stars = ifelse(pval < 0.001, "***",
                     ifelse(pval < 0.01, "**",
                            ifelse(pval < 0.05, "*",
                                   ifelse(pval < 0.1, "+ ", "n.s.")))),
      # Format p-value
      pval = ifelse(pval < 0.001, "< 0.001", sprintf("%.3f", pval))
    ) %>%
    filter(rowname != "(Intercept)") %>%
    mutate(term = term_labels,
           label = paste0(term, ": ", stars)
    )
  
  paste0(labs$label, collapse = "\n")
}

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
# dtd_n1 <- c("n = 90", '15', "79", '126')

diapauseData3 %>% group_by(treatment, core_edge) %>% summarise(
  mean = mean(time),
  sd = sd(time),
  n = n()
)


## model with core_Edge, pop as rand
## the poisson model is overdispersed, so going with a negative binomial model,
## which fixes the overdispersion
modDiapausing_pois <- glmmTMB(time ~ core_edge * treatment + (1| population), data = diapauseDataDiapausing,
                         family = poisson())
simulateResiduals(modDiapausing_pois) %>% plot()
testDispersion(modDiapausing_pois)

# Model - Diapausing, core_edge, pop as rand
# Negative binomial model to account for overdispersion
modDiapausing <- glmmTMB(time ~ core_edge * treatment + (1| population), data = diapauseDataDiapausing,
                         family = nbinom2())
simulateResiduals(modDiapausing) %>% plot()
testDispersion(modDiapausing)
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
# This doesn't totally fix overdispersion, but it is much better
modDiapausingpop <- glmmTMB(time ~  treatment * population, data = diapauseDataDiapausing,
                            family = nbinom2)
simulateResiduals(modDiapausingpop) %>% plot()
testDispersion(modDiapausingpop)
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

meansDiapausingPops_df <- meansDiapausingPops_df %>%
  left_join(., dtd_n2 %>% select(treatment, population, n, n2)) %>%
  mutate(letter = case_when(
    population == "Bl" ~ "H",
    population == "Lc" ~ "E",
    population == "Lj" ~ "G",
    population == "Wi" ~ "F",
    population == "De" ~ "C",
    population == "Hu" ~ "B",
    population == "Lo" ~ "A",
    population == "Pu" ~ "D"
  ))


## Plot of core/edge and population means


# Plot for Diapausers only with new format
diapausers_reaction_plot_new <- ggplot(data = meansDiapausing_df, 
                                       aes(x = treatment, y = response, color = core_edge)) +
  # add data points
  geom_point(data = diapauseDataDiapausing, aes(x = treatment, y = time, color = core_edge),
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.2), 
             alpha = 0.1) +
  #population lines and points
  geom_line(data = meansDiapausingPops_df, aes(group = population, color = core_edge),
            position = position_dodge(width = 0.4), alpha = 0.4) +
  geom_pointrange(data = meansDiapausingPops_df, 
                  aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL, 
                      group = population, color = core_edge),
                  position = position_dodge(width = 0.4), shape = 15, alpha = 0.6) +
  # core/edge lines and points
  geom_line(aes(group = core_edge, color = core_edge), 
            position = position_dodge(width = 0.3), linewidth = 0.8) +
  geom_pointrange(aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL),
                  position = position_dodge(width = 0.3), size = .7) +
  #sample size
  geom_text(aes(label = dtd_n1, y = -1.5), position = position_dodge(width = 0.5)) +
  # p-value labels
  annotate("text", label = plot_pval_labs(modDiapausing), 
           x = 2.5, y = 40, size = 3, hjust = 1, vjust = 1) +
  #formatting
  labs(x = "Daylength Treatment", y = "Days to Diapause", 
       # title = "Diapausing Individuals Only",
       color = "Origin") +
  scale_x_discrete(labels = c('long' = "Long\n(North)", 'short' = "Short\n(South)")) +
  scale_y_continuous(breaks = seq(0, 40, by = 10)) +
  scale_color_manual(values = c('core' = '#009E72', 'edge' = "#D55E00"), 
                     labels = c('edge' = "Edge", 'core' = "Core")) +
  guides(color = guide_legend(override.aes = list(shape = 19, linetype = 0))) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 13))
diapausers_reaction_plot_new


# Plot of just population differences 
diapausers_pop_plot_new <- ggplot(data = meansDiapausingPops_df) +
  ggpattern::geom_col_pattern(aes(x = treatment, y = response, 
                                  fill = letter, group = letter,
                                  pattern = core_edge),
                              pattern_density = 0.3,
                              pattern_fill = 'white',
                              position = position_dodge2(width = 1),
                              color = 'black',
                              alpha = 0.8) +
  ggpattern::scale_pattern_manual(values = c("none", "stripe"),
                                  labels = c('long' = "Long", 'short' = "Short")) + # Define pattern types
  geom_errorbar(data = meansDiapausingPops_df,
                aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL,
                    group = population), width = 0.9,
                position = position_dodge2(width = 2),
                color = 'black') +
  #sample size
  geom_text(data = meansDiapausingPops_df,
            aes(label = n, x = treatment, y = -1.5, group = letter),
            position = position_dodge2(width = 0.9),
            size = 2.9) +
  #formatting
  labs(x = "Daylength Treatment", y = "Days to Diapause", 
       # title = "Diapausing Individuals Only",
       fill = "Collection\nLocation", pattern = "Origin") +
  # scale_fill_manual(values = c('core' = '#009E72', 'edge' = "#D55E00"), 
  #                   labels = c('edge' = "Edge", 'core' = "Core")) +
  scale_fill_manual(values = colorRampPalette(colors = c("#009E72", "white", "#D55E00"))(8)) +
  scale_x_discrete(labels = c('long' = "Long\n(North)", 'short' = "Short\n(South)")) +
  guides(fill = guide_legend(override.aes = list(pattern = "none")),
         pattern = guide_legend(override.aes = list(fill = 'grey80'))) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 13))
diapausers_pop_plot_new

ggarrange(diapausers_reaction_plot_new, diapausers_pop_plot_new,
          common.legend = FALSE,
          legend = 'right', widths = c(.7, 1), labels = "AUTO") %>%
  annotate_figure(top = text_grob("Diapausing Individuals Only",
                                  size = 15))

ggsave("plots/S1_days_only_diapausers_new.jpg", 
       dpi = 400, width = 12, height = 6)


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

diapauseData3 %>%
  summarise(mean = mean(time),
            variance = var(time))

## model with core_Edge, pop as rand
# Poisson model is over dispersed
modAllSamp_pois <- glmmTMB(time ~ core_edge * treatment + (1| population), data = diapauseData3,
                      family = poisson())
simulateResiduals(modAllSamp) %>% plot()
testDispersion(modAllSamp)



# negative binomial model to fix over dispersion
modAllSamp <- glmmTMB(time ~ core_edge * treatment + (1| population), data = diapauseData3,
                      family = nbinom1())
simulateResiduals(modAllSamp) %>% plot()
testDispersion(modAllSamp)
summary(modAllSamp)
Anova(modAllSamp, type = 3)
emmeans(modAllSamp, pairwise ~ core_edge|treatment, type = 'response')
emmeans(modAllSamp, pairwise ~ treatment|core_edge, type = 'response')


emoutmodAllSamp <- emmeans(modAllSamp, pairwise ~ core_edge*treatment, type = "response")
emoutmodAllSamp_df <- as.data.frame(emoutmodAllSamp$emmeans)
emoutmodAllSamp_df$core_edge <- factor(emoutmodAllSamp_df$core_edge, levels = c('core', 'edge'))


# p-values for plotting
plot_pval_labs(modAllSamp)


# Population analyses - all samples
modAllSamppop <- glmmTMB(time ~  treatment * population, data = diapauseData3,
                         family = nbinom1())
simulateResiduals(modAllSamppop) %>% plot()
testDispersion(modAllSamppop)
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
emoutAllSampPops_df <- emoutAllSampPops_df %>% 
  left_join(., summ_dtd_all_pop %>% 
              select(treatment, population, n, n2))


## Combining population and core edge plots in reaction norm style


# New version of prop in diapause plot
days_plotnew <- ggplot(data = emoutmodAllSamp_df, aes(x = treatment, y = response, color = core_edge)) +
  # #box for non-diapausers
  annotate('rect', xmin = c(0.7), xmax = c(2.3), ymin = 41.7, ymax = 44.5, 
           fill = 'grey95', color = NA) +
  # non-diapauser label
  annotate("label", x = c(1.55), y = 43.3, label = "non-diapausers", 
           fill = NA, label.size = NA, size = 3.1) +
  # p-value labels
  annotate("text", label = plot_pval_labs(modAllSamp), 
           x = 0.5, y = 50, size = 3, hjust = 0, vjust = 1) +
  # add data points
  geom_point(data = diapauseData3, aes(x = treatment, y = time, color = core_edge),
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.2), 
             alpha = 0.1) +
  #population lines and points
  geom_line(data = emoutAllSampPops_df, aes(y = response, group = population, color = core_edge),
            position = position_dodge(width = 0.4), alpha = 0.6) +
  geom_pointrange(data = emoutAllSampPops_df, 
                  aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL, 
                      group = population, color = core_edge),
                  position = position_dodge(width = 0.4), shape = 15, alpha = 0.6) +
  # core/edge lines and points
  geom_line(aes(group = core_edge, color = core_edge), 
            position = position_dodge(width = 0.3)) +
  geom_pointrange(aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL),
                  position = position_dodge(width = 0.3), size = .7) +
  
  #sample size
  geom_text(data = summ_dtd_all, aes(label = n2, y = -1.5), position = position_dodge(width = 0.35)) +
    #formatting
  labs(x = "Daylength Treatment", y = "Days to Diapause", color = "Collection \nLocation") +
  scale_x_discrete(labels = c('long' = "Long\n(North)", 'short' = "Short\n(South)")) +
  scale_color_manual(values = c('core' = '#009E72', 'edge' = "#D55E00"), labels = c('edge' = "Edge", 'core' = "Core")) +
  scale_y_continuous(breaks = seq(0, 40, by = 10)) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())
days_plotnew



# Plot of just population differences 
dtd_all_pop_new <- ggplot(data = emoutAllSampPops_df) +
  ggpattern::geom_col_pattern(aes(x = treatment, y = response,
                                  fill = letter, group = letter,
                                  pattern = core_edge),
                              pattern_density = 0.3,
                              pattern_fill = 'white',
                              position = position_dodge2(width = 1),
                              color = 'black',
                              alpha = 0.8) +
  ggpattern::scale_pattern_manual(values = c("none", "stripe"),
                                  labels = c('long' = "Long", 'short' = "Short")) + # Define pattern types
  geom_errorbar(data = emoutAllSampPops_df,
                aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL,
                    group = population), width = 0.9,
                position = position_dodge2(width = 2),
                color = 'black') +
  #sample size
  geom_text(data = emoutAllSampPops_df,
            aes(label = n, x = treatment, y = -1.5, group = letter),
            position = position_dodge2(width = 0.9),
            size = 2.9) +
  #formatting
  labs(x = "Daylength Treatment", y = "Days to Diapause",
       title = "Days to Diapause - including non-diapausers",
       fill = "Collection\nLocation", pattern = "Origin") +
  # scale_fill_manual(values = c('core' = '#009E72', 'edge' = "#D55E00"),
  #                   labels = c('edge' = "Edge", 'core' = "Core")) +
  scale_fill_manual(values = colorRampPalette(colors = c("#009E72", "white", "#D55E00"))(8)) +
  scale_x_discrete(labels = c('long' = "Long\n(North)", 'short' = "Short\n(South)")) +
  guides(fill = guide_legend(override.aes = list(pattern = "none")),
         pattern = guide_legend(override.aes = list(fill = 'grey80'))) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 13))


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
survreg_means_df <- survreg_means_df %>% 
  mutate(
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
survreg_means_df <- survreg_means_df %>% 
  left_join(., summ_dtd_all_pop %>% 
              select(treatment, population, n, n2))

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


# # New version of prop in diapause plot
surv_plot_new <- ggplot(data = survreg_means_CE_df, aes(x = treatment, y = emmean, color = core_edge)) +
  # #box for non-diapausers
  annotate('rect', xmin = c(0.7), xmax = c(2.3), ymin = 41.7, ymax = 44.5,
           fill = 'grey95', color = NA) +
  # non-diapauser label
  annotate("label", x = c(1.55), y = 43.3, label = "non-diapausers",
           fill = NA, label.size = NA, size = 3.1) +
  # add data points
  geom_point(data = diapauseData3, aes(x = treatment, y = time, color = core_edge),
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.2),
             alpha = 0.1) +
  #population lines and points
  geom_line(data = survreg_means_df, aes(group = population, color = core_edge),
            position = position_dodge(width = 0.4), alpha = 0.6) +
  geom_pointrange(data = survreg_means_df,
                  aes(x = treatment, y = emmean, ymin = lower.CL, ymax = upper.CL,
                      group = population, color = core_edge),
                  position = position_dodge(width = 0.4), shape = 15, alpha = 0.6) +
  # core/edge lines and points
  geom_line(aes(group = core_edge, color = core_edge),
            position = position_dodge(width = 0.3),
            linewidth = 0.8) +
  geom_pointrange(aes(x = treatment, y = emmean, ymin = lower.CL, ymax = upper.CL),
                  position = position_dodge(width = 0.3), size = .7) +
  #sample size
  geom_text(data = summ_dtd_all, aes(label = n2, y = -5), position = position_dodge(width = 0.35)) +
  #formatting
  labs(x = "Daylength Treatment", 
       y = "Days to Diapause", 
       # title = "Days to Diapause - Survival Analysis",
       color = "Origin") +
  scale_x_discrete(labels = c('long' = "Long\n(North)", 'short' = "Short\n(South)")) +
  scale_color_manual(values = c('core' = '#009E72', 'edge' = "#D55E00"), 
                     labels = c('edge' = "Edge", 'core' = "Core")) +
  guides(color = guide_legend(override.aes = list(shape = 19, linetype = 0))) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 13))


 
# # Plot of just population differences 
surv_pop_plot_new <- ggplot(data = survreg_means_df) +
  ggpattern::geom_col_pattern(aes(x = treatment, y = emmean,
                                  fill = letter, group = letter,
                                  pattern = core_edge),
                              pattern_density = 0.3,
                              pattern_fill = 'white',
                              position = position_dodge2(width = 1),
                              color = 'black',
                              alpha = 0.8) +
  ggpattern::scale_pattern_manual(values = c("none", "stripe"),
                                  labels = c('long' = "Long", 'short' = "Short")) + # Define pattern types
  geom_errorbar(data = survreg_means_df,
                aes(x = treatment, ymin = lower.CL, ymax = upper.CL,
                    group = population), width = 0.9,
                position = position_dodge2(width = 2),
                color = 'black') +
  #sample size
  geom_text(data = survreg_means_df,
            aes(label = n, x = treatment, y = upper.CL + 3, group = letter),
            position = position_dodge2(width = 0.9),
            size = 2.9) +
  #formatting
  labs(x = "Daylength Treatment", y = "Days to Diapause",
       # title = "Diapausing Individuals Only",
       fill = "Collection\nLocation", pattern = "Origin") +
  # scale_fill_manual(values = c('core' = '#009E72', 'edge' = "#D55E00"),
  #                   labels = c('edge' = "Edge", 'core' = "Core")) +
  scale_fill_manual(values = colorRampPalette(colors = c("#009E72", "white", "#D55E00"))(8)) +
  scale_x_discrete(labels = c('long' = "Long\n(North)", 'short' = "Short\n(South)")) +
  guides(fill = guide_legend(override.aes = list(pattern = "none")),
         pattern = guide_legend(override.aes = list(fill = 'grey80'))) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 13))

ggarrange(surv_plot_new, surv_pop_plot_new, common.legend = FALSE,
          legend = 'right', labels = "AUTO", widths = c(0.7, 1)) %>%
  annotate_figure(top = text_grob("Days to Diapause - Survival Analysis",
                  size = 15))

ggsave("plots/S2_days_survival_new.png", dpi = 400, width = 12, height = 6)


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
logmodel1 <- glmmTMB(diapause ~ core_edge * treatment + (1|population), 
                     family = "binomial", data = diapauseData3)
simulateResiduals(logmodel1) %>% plot()
testDispersion(logmodel1)
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
logmodel_bayes <- arm::bayesglm(diapause ~ population * treatment, 
                                family = 'binomial', data = diapauseData3)
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
logisticbayes_means_pop_df <- logisticbayes_means_pop_df %>%
  left_join(., prop_sample_size ) %>%
  mutate(n2 = replace(n, n==26, 'n=26') )

plot_pval_labs(logmodel_bayes, 
               term_labels = c("Population", "Treatment", "Population * Treatment"))

## Combining population and core edge plots in reaction norm style

# New version of prop in diapause plot
prop_plotnew <- ggplot(data = logistic_means_ce_df, aes(x = treatment, y = prob)) +
  #population lines and points
  geom_line(data = logisticbayes_means_pop_df,
            aes(x = treatment, y = prob,
                group = population, color = core_edge),
            position = position_dodge(width = 0.4),
            alpha = 0.3) +
  geom_pointrange(data = logisticbayes_means_pop_df,
                  aes(x = treatment, y = prob, ymin = asymp.LCL, ymax = asymp.UCL,
                      group = population, color = core_edge),
                  position = position_dodge(width = 0.4),
                  alpha = 0.3, shape = 15) +
  annotate('text', label = plot_pval_labs(logmodel1),
           x = 0.5, y = 1.05, hjust = 0, vjust = 1,
           size = 3) +
  #core/edge lines and points
  geom_line(aes(color = core_edge, group = core_edge), 
            position = position_dodge(width = 0.3)) +
  geom_pointrange(aes(x = treatment, y = prob, ymin = asymp.LCL, ymax = asymp.UCL, 
                      color = core_edge),
                  position = position_dodge(width = 0.3), size = 0.7) +
  #sample size
  geom_text(data = summ_dtd_all, aes(label = n2, y = 0, color = core_edge), position = position_dodge(width = 0.4)) +
  #formatting
  labs(x = "Daylength Treatment", y = "Proportion in Diapause", color = "Origin") +
  scale_x_discrete(labels = c('long' = "Long\n(North)", 'short' = "Short\n(South)")) +
  scale_color_manual(values = c('core' = '#009E72', 'edge' = "#D55E00"), labels = c('edge' = "Edge", 'core' = "Core")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  guides(color = guide_legend(override.aes = list(shape = 19, linetype = 0))) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())
prop_plotnew


ggarrange(prop_plotnew, days_plotnew, nrow = 1, labels = "AUTO", common.legend = T, legend = 'bottom')

ggsave("plots/Fig2_prop_days_new.png", dpi = 400, width = 10, height = 6)


# Plot only populations - proportion
prop_pop_new <- ggplot(data = logisticbayes_means_pop_df) +
  ggpattern::geom_col_pattern(aes(x = treatment, y = prob,
                                  fill = letter, group = letter,
                                  pattern = core_edge),
                              pattern_density = 0.3,
                              pattern_fill = 'white',
                              position = position_dodge2(width = 1),
                              color = 'black',
                              alpha = 0.8) +
  ggpattern::scale_pattern_manual(values = c("none", "stripe"),
                                  labels = c('long' = "Long", 'short' = "Short")) + # Define pattern types
  geom_errorbar(data = logisticbayes_means_pop_df,
                aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL,
                    group = population), width = 0.9,
                position = position_dodge2(width = 2),
                color = 'black') +
  #sample size
  geom_text(data = logisticbayes_means_pop_df,
            aes(label = n, x = treatment, y = -0.05, group = letter),
            position = position_dodge2(width = 0.9),
            size = 2.9) +
  #formatting
  labs(x = "Daylength Treatment", y = "Proportion in Diapause",
       title = "Proportion in Diapause",
       fill = "Collection\nLocation", pattern = "Origin") +
  # scale_fill_manual(values = c('core' = '#009E72', 'edge' = "#D55E00"),
  #                   labels = c('edge' = "Edge", 'core' = "Core")) +
  scale_fill_manual(values = colorRampPalette(colors = c("#009E72", "white", "#D55E00"))(8)) +
  scale_x_discrete(labels = c('long' = "Long\n(North)", 'short' = "Short\n(South)")) +
  guides(fill = guide_legend(override.aes = list(pattern = "none")),
         pattern = guide_legend(override.aes = list(fill = 'grey80'))) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 13))

ggarrange(prop_pop_new, dtd_all_pop_new, common.legend = TRUE,
          legend = 'right', labels = "AUTO")

ggsave("plots/S2_pop_figs.png", 
       dpi = 400, width = 10, height = 5)

# Tables of results ------------------------------------------------------------

####
#' Format ANOVA Table for Publication
#'
#' Creates a publication-ready ANOVA table using tinytable from a fitted model
#'
#' @param model A fitted model object
#' @param term_names Optional character vector of custom term names. If NULL, uses rownames from ANOVA table
#' @param caption Optional caption for the table. If NULL, uses a default caption
#' @param notes Optional notes for the table. If NULL, uses default significance codes
#' @param anova_type Type of ANOVA (default = 3)
#'
#' @return A formatted tinytable object
#'
#' @examples
#' # With default term names
#' format_anova_table(modDiapausing)
#' 
#' # With custom term names
#' format_anova_table(modDiapausing, 
#'                    term_names = c("Intercept", "Origin", "Treatment", "Origin * Treatment"))
#'
format_anova_table <- function(model, 
                               term_names = NULL,
                               caption = NULL,
                               notes = "Significance codes: *** P < 0.001, ** P < 0.01, * P < 0.05, † P < 0.1",
                               anova_type = 3) {
  
  # Get ANOVA table
  anova_results <- Anova(model, type = anova_type)
  
  # Convert to data frame and prepare for formatting
  anova_df <- as.data.frame(anova_results)
  
  # Use custom term names if provided, otherwise use rownames
  if (!is.null(term_names)) {
    if (length(term_names) != nrow(anova_df)) {
      stop("Length of term_names must match number of rows in ANOVA table")
    }
    anova_df$Term <- term_names
  } else {
    anova_df$Term <- rownames(anova_df)
  }
  
  anova_df <- anova_df[, c("Term", "Chisq", "Df", "Pr(>Chisq)")]
  
  # Set default caption if not provided
  if (is.null(caption)) {
    caption <- paste0("Type ", anova_type, " Analysis of Deviance Table (Wald χ² tests)")
  }
  
  # Create publication-ready table
  tt(anova_df, 
     caption = caption,
     notes = notes) |>
    # Format column names
    setNames(c("Term", "χ²", "df", "P-value")) |>
    # Format numeric columns
    format_tt(digits = 3, num_fmt = "decimal") |>
    # Format p-values and add significance stars in one step
    format_tt(j = 4, fn = function(x) {
      # Determine stars
      stars <- ifelse(x < 0.001, "***",
                      ifelse(x < 0.01, "**",
                             ifelse(x < 0.05, "*",
                                    ifelse(x < 0.1, "†", ""))))
      # Format p-value
      p_formatted <- ifelse(x < 0.001, "< 0.001", sprintf("%.3f", x))
      # Combine with stars (only add space if stars exist)
      ifelse(stars == "", p_formatted, paste0(p_formatted, " ", stars))
    }) |>
    # Style the table
    style_tt(
      i = 0,  # Header row
      bold = TRUE,
      line = "b",
      line_width = 0.2
    ) |>
    style_tt(
      i = nrow(anova_df),  # Last row
      line = "b",
      line_width = 0.2
    )
}
 
# create tables of results for two neg. binom modles & logistic model
diaptable <- format_anova_table(modDiapausing,
                   term_names = c("Intercept", "Origin", "Treatment", "Origin * Treatment"))
save_tt(diaptable, output = "plots/diap_table.png")
alltable <- format_anova_table(modAllSamp,
                   term_names = c("Intercept", "Origin", "Treatment", "Origin * Treatment"))
save_tt(alltable, "plots/alltable.png")
proptable <- format_anova_table(logmodel1,
                   term_names = c("Intercept", "Origin", "Treatment", "Origin * Treatment"))
save_tt(proptable, "plots/proptable.png")


# create a table using tt of the results from logisticbayes_means_pop_df 
logisticbayes_means_pop_df %>%
  mutate(
    # Format probability and confidence intervals
    `Probability (95% CI)` = sprintf("%.3f (%.3f–%.3f)", 
                                     prob, asymp.LCL, asymp.UCL),
    # # Format SE
    # SE = sprintf("%.3f", SE),
    # Clean up treatment names
    Treatment = ifelse(treatment == "long", "Long", "Short"),
    # # Rename population for clarity
    # Population = population,
    # Keep core/edge and letter groupings
    Origin = ifelse(core_edge == 'core', "Core", "Edge"),
    Population = letter
  ) %>%
  select(Origin, Treatment, Population, `Probability (95% CI)`) %>%
  arrange(Origin, Treatment, Population) %>%
  # Create the tinytable
  tt()

emoutAllSampPops_df %>%
  mutate(
    # Format probability and confidence intervals
    `Days to diapause (95% CI)` = sprintf("%.3f (%.3f–%.3f)", 
                                     rate, asymp.LCL, asymp.UCL),
    # # Format SE
    # SE = sprintf("%.3f", SE),
    # Clean up treatment names
    Treatment = ifelse(treatment == "long", "Long", "Short"),
    # # Rename population for clarity
    # Population = population,
    # Keep core/edge and letter groupings
    Origin = ifelse(core_edge == 'core', "Core", "Edge"),
    Population = letter
  ) %>%
  select(Origin, Treatment, Population, `Days to diapause (95% CI)`) %>%
  arrange(Origin, Treatment, Population) %>%
  # Create the tinytable
  tt()


##Correlation between prop & days ---------------------------------------------


##Using population means

logisticbayes_means_pop_df

emoutAllSampPops_df

prop_days_pop <- full_join(logisticbayes_means_pop_df,emoutAllSampPops_df, by = c('treatment', 'population'))
corr_diapause <- cor.test(prop_days_pop$prob, prop_days_pop$rate, method = c("pearson"))

ggplot(data = prop_days_pop)+
  geom_smooth(aes(x = prob, y = rate), method = lm, color = 'black') +
  geom_point(aes(x = prob, y = rate, fill = core_edge.x, shape = treatment),
              size = 3) +
  # geom_text(aes(label = letter, x = prob+.02, y = rate + 1), size = 4) +
  labs(y = "Days until Diapause", x = "Proportion in Diapause", fill = 'Origin',
       shape = 'Treatment') +
  annotate("text", y = 40, x = 0.75, 
           label = paste0("italic(R) ^ 2 == ", round(corr_diapause$estimate,2)), parse = TRUE, size = 6) +
  annotate("text", y = 37, x = 0.75, 
           label = "italic(P) < 0.001", parse = TRUE, size = 5) +
  # annotate("text", y = -2, x = site_char$latitude[1:4], label = site_char$letter[1:4], size = 4) +
  scale_shape_manual(values = c(21, 22), 
                     labels = c('long' = "Long (North)", 'short' = "Short (South)")) +
  scale_fill_manual(values = c('core' = '#009E72', 'edge' = "#D55E00"), 
                    labels = c('edge' = "Edge", 'core' = "Core")) +
  guides(fill = guide_legend('Origin', override.aes = list(shape = 21))) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

ggsave("plots/S3_corr_days_prop.png", dpi = 400, width = 7, height = 4.5)


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



### Version 2 -------------------------------

# PCA on elevation, latitude, degree days
# Perform PCA on numeric variables (excluding population and letter)
pca_vars <- site_char[, c("latitude", "elevation", "degree_days", "first_frost")]
pca_result <- prcomp(pca_vars, scale. = TRUE, center = TRUE)

# Print PCA summary
summary(pca_result)
pca_result$rotation

# Extract PC1 scores for each population
pc1_scores <- data.frame(
  population = site_char$population,
  letter = site_char$letter,
  PC1 = pca_result$x[, 1]
)


# Create PCA plot
pca_plot_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  population = site_char$population,
  core_edge = factor(site_char$core_edge, labels = c("core", "edge")),
  letter = site_char$letter
)

# Calculate variance explained
var_explained <- round(100 * summary(pca_result)$importance[2, ], 1)

# Prepare biplot data with importance-weighted arrows
loadings <- pca_result$rotation[, 1:2]

# Calculate the importance of each variable (contribution to variance)
# This is based on the squared loadings weighted by variance explained
importance <- sqrt(loadings[, 1]^2 * var_explained[1] + 
                     loadings[, 2]^2 * var_explained[2])

# Scale loadings by their importance and add a visibility scaling factor
scale_factor <- .25
loadings_scaled <- loadings * importance * scale_factor

loading_data <- data.frame(
  variable = rownames(loadings),
  PC1 = loadings_scaled[, 1],
  PC2 = loadings_scaled[, 2],
  importance = importance
) %>%
  mutate(labels = gsub(variable, pattern = "_", replacement = " "))


# PCA biplot with loadings and populations
ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = core_edge, label = letter)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_text(hjust = -1, vjust = 0.5, size = 3.5, show.legend = FALSE) +
  # Add variable loading arrows
  geom_segment(data = loading_data, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "black", size = 0.8, inherit.aes = FALSE) +
  geom_text(data = loading_data,
            aes(x = PC1, y = PC2, label = labels),
            color = "black", size = 3.5,
            hjust = ifelse(loading_data$PC1 > 0, -0.2, 1.2),
            vjust = ifelse(loading_data$PC2 > 0, -0.5, 1.5),
            inherit.aes = FALSE) +
  labs(
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)"),
    color = "Location"
  ) +
  lims(x = c(-3.8, 2.3)) +
  scale_color_manual(values = c('core' = '#009E72', 'edge' = "#D55E00")) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        legend.position = 'bottom',
        # plot.title = element_text(hjust = 0.5, face = "bold"),
        ) 

ggsave("plots/PCA_biplot.png", dpi = 400, width = 6, height = 5)





# Merge PC1 scores with diapause data
diapauseData_pc1 <- diapauseData3 %>% 
  mutate(lat_sc = scale(latitude)[,1],
         ele_sc = scale(elevation)[,1],
         deg_sc = scale(degree_days)[,1]) %>%
  left_join(pc1_scores, by = "population")




# models with PC1
library(lme4)
library(lmerTest)
mod_PC1_all <- lmer(time ~ PC1 * treatment + (1|population), 
                       data = diapauseData_pc1)
summary(mod_PC1_all)
anova(mod_PC1_all)
MuMIn::r.squaredGLMM(mod_PC1_all)
emtrends(mod_PC1_all, ~ treatment, var = "PC1")

# # quadratic model
# mod_PC1_quad <- lmer(time ~ PC1 * treatment * I(PC1^2) + (1|population), 
#                     data = diapauseData_pc1)
# summary(mod_PC1_quad)
# AIC(mod_PC1_all, mod_PC1_quad)


# Make predictions 

pc1_range <- seq(min(diapauseData_pc1$PC1), 
                 max(diapauseData_pc1$PC1), 
                 length.out = 100)

# Get unique treatment levels
treatment_levels <- unique(diapauseData_pc1$treatment)

# Create prediction grid for both treatments
pred_data <- expand.grid(
  PC1 = pc1_range,
  treatment = treatment_levels
)

# Make predictions at population level (marginal predictions)
# This averages over random effects
pred_data$predicted_time <- predict(mod_PC1_all, 
                                    newdata = pred_data, 
                                    re.form = NA)  # Exclude random effects

# # Calculate confidence intervals using parametric bootstrap or standard errors
# # Using the fixed effects only
mm <- model.matrix(~ PC1 * treatment, data = pred_data)
# pred_data$predicted_time_fe <- mm %*% fixef(mod_PC1_all)

# Get variance-covariance matrix for standard errors
vcov_mat <- vcov(mod_PC1_all)
pred_data$se <- sqrt(diag(mm %*% vcov_mat %*% t(mm)))

# Calculate 95% confidence intervals
pred_data$lower <- pred_data$predicted_time - 1.96 * pred_data$se
pred_data$upper <- pred_data$predicted_time + 1.96 * pred_data$se



labs = c("Long Daylength", "Short Daylength")
names(labs) = c("long", "short")

trend_labs <- emtrends(mod_PC1_all, ~ treatment, var = "PC1") %>% 
  as.data.frame() %>%
  mutate(label = paste0("Slope = ", round(PC1.trend, 2), 
                        "\n(", round(lower.CL, 2), ", ",
                        round(upper.CL, 2), ")"))

# Plot
ggplot(data = diapauseData_pc1, 
       aes(x = PC1, y = time, group = treatment)) +
  geom_jitter(aes(fill = letter), width = .1, height = 0.2, shape = 21, size = 3, alpha = 0.5) +
  geom_ribbon(data = pred_data, aes(x = PC1, ymin = lower, ymax = upper), 
              inherit.aes = FALSE, alpha = 0.2) +
  geom_line(data = pred_data, aes(x = PC1, y = predicted_time), 
            color = 'black', size = 1) +
  geom_text(data = trend_labs, 
            aes(group = treatment, label = label), 
            x = -3.5, y = -5, hjust = 0, size = 3.5) +
  # geom_smooth(method = lm, color = 'blue') +
  facet_grid(~ treatment, labeller = labeller(treatment = labs)) +
  labs(y = "Days to Diapause", x = "Environment (PC1)", 
       fill = "Collection\nLocation") +
  # scale_fill_brewer(palette = "PiYG", direction = -1) +
  scale_y_continuous(breaks = seq(0, 50, by = 10)) +
  scale_fill_manual(values = colorRampPalette(colors = c("#009E72", "white", "#D55E00"))(8)) +
  # scale_fill_manual(values = c('core' = '#414487FF', 'edge' = "#7AD151FF"), labels = c('edge' = "South", 'core' = "North")) +
  # guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

ggsave("plots/env_correlations_new.png", dpi = 400, width = 10, height = 5)

