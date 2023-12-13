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

# Data load and clean up ----
#load diapause data, delete missing/messed up trials, create time column, make trt & pop factors
diapauseData <- read_excel("diapause.xlsx")
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

  
#load fecundity/weight data
weightData <- read_xlsx("G:/My Drive/GRAD SCHOOL/RESEARCH/D. carinulata/range expansion dynamics experiments/eggcount.xlsx")

#Combine data sets
diapauseData3 <- left_join(diapauseData2, 
          weightData %>% select(beetle_ID, eggs, weight), 
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
  mutate(distFromCore = 39.144 - latitude) %>%
  mutate(distFromEdge = 34.422 - latitude)
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


# PCA of climate variables ----
site_char <- data.frame(
  population = c("Bl", "Lc", "Lj", "Wi", "De", "Hu", "Lo", "Pu"),
  latitude = c(33.912,34.593,34.342,34.422,39.144,40.063,44.856,38.268),
  elevation = c(96.9, 1669.56, 1433, 1204.74,1386.43, 1189.62, 1115.35, 1448.19),
  degree_days = c(4524, 2043, 2547, 3516, 1767, 1841, 1438, 1990),
  core_edge = c(rep("1", 4), rep("0", 4))
  )
site_char$core_edge <- as.numeric(site_char$core_edge)
str(site_char)

biplot(prcomp(site_char[,2:5], scale = TRUE), scale = 0)


# Objects for plots, etc.

edge_col = viridis(3)[1]
core_col = viridis(3)[2]

# SURVIVAL ANALYSES ----
## Create survival object
survObj <- Surv(diapauseData3$time, diapauseData3$diapause)

# Kaplan-Meier model - population ----
model1 <- survfit(survObj ~ population + treatment, data = diapauseData3)
model1
summary(model1)
model1_summary <- surv_summary(model1, data = diapauseData3) %>%
  select(time, surv, population, treatment)

time_zero = data.frame(time = rep(0, 16),
                       surv = rep(1, 16),
                       population = rep(c("Bl", "Lc", "Lj", "Wi", "De", "Hu", "Lo", "Pu"), each = 2),
                       treatment = rep(c("long", "short"),8))

model1_summary2 <- rbind(model1_summary, time_zero) %>%
  mutate(core_edge = case_when(
  population == "Bl" | population == "Lc" | population == "Lj" | population == "Wi" ~ "edge",
  population == "De" | population == "Hu" | population == "Lo" | population == "Pu" ~ "core")) %>%
  mutate(invSurv = 1-surv)

## Plots of survival curves----
ggplot(data = model1_summary2, aes(x = time, y = invSurv , color = core_edge, fill = population)) +
  geom_step(size = .5, position = position_dodge(width = 0.8)) +
  facet_grid(~ treatment) +
  scale_color_manual(values = c(core_col, edge_col)) +
  # geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.1, linetype = 0, color="grey") +
  labs(y = "Probability of diapause", x = "Time (days)") +
  theme_bw() +
  theme(legend.position = "none")
# ggsave(filename = "adapt both.svg", path = "G:/My Drive/GRAD SCHOOL/CONFERENCES/Entomology2021",
       # plot = last_plot(), device = "svg", width = 5.36, height = 2.74, units = "in")

#same plot, only core
ggplot(data = model1_summary2 %>% filter(core_edge == "core"), aes(x = time, y = invSurv , color = core_edge, fill = population)) +
  geom_step(size = .5, position = position_dodge(width = 0.8)) +
  facet_grid(~ treatment) +
  scale_color_manual(values = c(core_col)) +
  # geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.1, linetype = 0, color="grey") +
  labs(y = "Probability of diapause", x = "Time (days)") +
  theme_bw() +
  theme(legend.position = "none")
# ggsave(filename = "adapt core.svg", path = "G:/My Drive/GRAD SCHOOL/CONFERENCES/Entomology2021",
       # plot = last_plot(), device = "svg", width = 5.36, height = 2.74, units = "in")

#same plot, only edge
ggplot(data = model1_summary2 %>% filter(core_edge == "edge"), aes(x = time, y = invSurv , color = core_edge, fill = population)) +
  geom_step(size = .5, position = position_dodge(width = 0.8)) +
  facet_grid(~ treatment) +
  scale_color_manual(values = c(edge_col)) +
  # geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.1, linetype = 0, color="grey") +
  labs(y = "Probability of diapause", x = "Time (days)") +
  theme_bw() +
  theme(legend.position = "none")
# ggsave(filename = "adapt edge.svg", path = "G:/My Drive/GRAD SCHOOL/CONFERENCES/Entomology2021",
       # plot = last_plot(), device = "svg", width = 5.36, height = 2.74, units = "in")

#default survival plot
ggsurvplot(model1, data = diapauseData3, 
           conf.int = F,
           surv.median.line = "hv",
           facet.by = "treatment",
           legend.title = "population",
           # legend.labs = c("north", "south"),
           ylab = "Probability not in diapause",
           xlab = "Time (days)",
           censor.shape = 124,
           censor = F,
           ggtheme = theme_bw())
           # legend = "right")

# Survival Regression Models ---- 
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
  unite("treatment_core_edge", c(treatment,core_edge), remove = F) %>%
  add_row(treatment_core_edge = 'long_x', treatment = 'long', population = 'XX', emmean = NA, SE = NA, df = NA, lower.CL = NA, upper.CL = NA, core_edge = 'core')

## Core and Edge
model1b_2 <- survreg(survObj ~ core_edge * treatment + population, dist ="gaussian", data = diapauseData3)
model1b_2
summary(model1b_2)
emmeans(model1b_2, pairwise ~ treatment | core_edge, type = "response")
survreg_means_CE <- emmeans(model1b_2, pairwise ~ treatment | core_edge, type = "response")
survreg_means_CE_df <- as.data.frame(survreg_means_CE$emmeans) %>%
  unite("treatment_core_edge", c(treatment,core_edge), remove = F)

# Plot for local adaptation with lines
grayscale_gradient <- c("gray25", "gray90", "gray75", "gray50", "gray10")
color_gradient <- c("red", "green", "yellow", "red", "darkred")

g <- make_gradient(
  deg = 270, n = 500, cols = grayscale_gradient
)

ggplot(data = survreg_means_CE_df, aes(x = treatment, y = emmean)) +
  # annotation_custom(
  #   grob = g, xmin = 2.4, xmax = Inf, ymin = 0, ymax = Inf) +
  geom_pointrange(aes(ymax = emmean + SE, ymin = emmean - SE, color = core_edge), 
                  size = 1.2, fatten = .85, stroke = 3, lty = 'solid') +
  labs(y = "Days until Diapause", x = "Fall Environment", color = "Collection \nLocation") +
  scale_color_manual(values = c(edge_col, core_col), labels = c('core' = 'North', 'edge' = 'South')) +
  scale_x_discrete(labels = c('long' = "North", 'short' = "South"))+
  # scale_y_continuous(limits = c(0,1.3)) +
  geom_line(data = survreg_means_df, aes(color = core_edge, group = population), 
            position = position_dodge(width = 0.4), show.legend = F) +
  geom_pointrange(data = survreg_means_df,
                  aes(x = treatment, y = emmean, ymax = emmean + SE, ymin = emmean - SE, color = core_edge, group = population), 
                  shape = 17, size = 0.9, alpha = 0.5, show.legend = T,
                  position = position_dodge(width = 0.4) ) +
  theme_classic2() +
  theme(text = element_text(size = 25), axis.ticks.x = element_blank())


############################
#Survival Plot for paper - revision -----
survreg_means_df$population <- factor(survreg_means_df$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))


ggplot(data = survreg_means_CE_df, aes(x = treatment, y = emmean)) +
  # annotation_custom(
  #   grob = g, xmin = 2.4, xmax = Inf, ymin = 0, ymax = Inf) +
  geom_pointrange(aes(ymax = emmean + SE, ymin = emmean - SE, color = core_edge), 
                  size = .51, stroke = 1, shape = 15, lty = 'solid') +
  labs(y = "Days until Diapause", x = "Fall Environment", color = "Collection \nLocation") +
  scale_color_manual(values = c(edge_col, core_col), labels = c('core' = 'North', 'edge' = 'South')) +
  scale_x_discrete(labels = c('long' = "North", 'short' = "South"))+
  # scale_y_continuous(limits = c(0,1.3)) +
  geom_line(data = survreg_means_df, aes(color = core_edge, group = population), position = position_jitter(width = 0.2, seed = 5), show.legend = F) +
  geom_pointrange(data = survreg_means_df,
                  aes(x = treatment, y = emmean, ymax = emmean + SE, ymin = emmean - SE, color = core_edge),
                  position = position_jitter(width = 0.2, seed = 5), shape = 16, size = 1, alpha = 0.5, show.legend = T) +
  theme_classic2() +
  theme(text = element_text(size = 25), axis.ticks.x = element_blank())

##Plot for paper - revision in style of other plots in paper
surv_coreedge_plot <- ggplot(data = survreg_means_CE_df, aes(x = treatment, y = emmean, fill = core_edge)) +
  # annotate('rect', xmin = c(0.8,1.8), xmax = c(1.2, 2.2), ymin = 42, ymax = 44, alpha = 0.2, color = "black") +
  # geom_point(data = diapauseData3, aes(x = treatment, y = time_diap, fill = core_edge),
  #            position = position_jitterdodge(dodge.width = 0.35, jitter.width = 0.15),
  #            alpha = 0.25) +
  geom_line(aes(group = core_edge), position = position_dodge(width = 0.35), size = 0.73) +
  geom_pointrange(aes(x = treatment, y = emmean, ymin = lower.CL, ymax = upper.CL),
                  position = position_dodge(width = 0.35), shape = 21, size = 0.8) +
  # geom_text(data = summ_dtd_all, aes(label = n2, y = -1.5), position = position_dodge(width = 0.35)) +
  # annotate("label", x = c(1,2), y = 45.5, label = "Non-Diapausers") +
  labs(x = "Diapause-Inducing Daylength", y = "Days to Diapause (95% CI)", fill = "Collection \nLocation") +
  scale_x_discrete(labels = c('long' = "North", 'short' = "South")) +
  scale_fill_manual(values = c('core' = 'white', 'edge' = 'black'), labels = c('edge' = "South", 'core' = "North")) +
  theme_bw(base_size = 15)
surv_coreedge_plot

surv_pop_plot <- ggplot(data = survreg_means_df %>% filter(treatment_core_edge != 'long_x')) +
  # annotate('rect', xmin = 0.75, xmax = 8.25, ymin = 42, ymax = 44, alpha = 0.2, color = 'black') +
  # geom_point(data = diapauseData3, aes(x = population, y = time_diap),
  #            position = position_jitter(width = 0.15), alpha = 0.25) +
  geom_pointrange(aes(x = population, y = emmean, ymin = lower.CL, ymax = upper.CL, fill = core_edge),
                  color = 'black', shape = 21, size = 0.8) +
  facet_grid( ~ treatment, labeller = labeller(treatment = c('long' = "Northern Daylength", 'short' = "Southern Daylength"))) +
  # geom_text(data = summ_dtd_all_pop, aes(label = n2, x = population, y = -1)) +
  # annotate("label", x = 4.5, y = 46, label = "Non-Diapausers") +
  labs(x = "Collection Location (N to S)", y = "Days to Diapause (95% CI)", fill = "Collection \nLocation") +
  scale_x_discrete(labels = c('Lo' = 'A', 'Hu' = 'B', 'De' = 'C', 'Pu' = 'D', 'Lc' = 'E', 'Wi' = 'F', 'Lj' = 'G', 'Bl' = 'H')) +
  scale_fill_manual(values = c('core' = 'white', 'edge' = 'black'), labels = c('edge' = "South", 'core' = "North")) +
  theme_bw(base_size = 15) +
  theme(strip.background = element_rect(fill = 'white'))
surv_pop_plot

ggarrange(surv_coreedge_plot, surv_pop_plot, 
          common.legend = T, 
          legend = 'right', 
          widths = c(1.3,2),
          labels = 'AUTO')
##export: width = 1100, height = 575 (11.46in x 5.99in)


##Plot for EntSoc presentation 2022
ggplot(data = survreg_means_CE_df, aes(x = treatment, y = emmean)) +
  geom_pointrange(aes(ymax = emmean + SE, ymin = emmean - SE, color = core_edge), 
                  size = 1.2, stroke = 2, lty = 'solid', show.legend = F) +
  labs(y = "Days until Diapause", x = "Fall Environment", color = "Collection \nLocation") +
  scale_color_manual(values = c(edge_col, core_col), labels = c('core' = 'North', 'edge' = 'South')) +
  scale_x_discrete(labels = c('long' = "North", 'short' = "South"))+
  # scale_y_continuous(limits = c(0,1.3)) +
  geom_line(data = survreg_means_df, aes(color = core_edge, group = population), position = position_jitter(width = 0.2, seed = 5), show.legend = F) +
  geom_pointrange(data = survreg_means_df,
                  aes(x = treatment, y = emmean, ymax = emmean + SE, ymin = emmean - SE, color = core_edge),
                  position = position_jitter(width = 0.2, seed = 5), size = 0.9, alpha = 0.5, show.legend = T) +
  theme_classic2() +
  theme(text = element_text(size = 25), axis.ticks.x = element_blank())
#without lines
ggplot(data = survreg_means_CE_df, aes(x = treatment, y = emmean)) +
  geom_pointrange(aes(ymax = emmean + SE, ymin = emmean - SE, color = core_edge), 
                  size = 1.2, stroke = 2, lty = 'solid', show.legend = F) +
  labs(y = "Days until Diapause", x = "Fall Environment", color = "Collection \nLocation") +
  scale_color_manual(values = c('white', core_col), labels = c('core' = 'North', 'edge' = 'South')) +
  scale_x_discrete(labels = c('long' = "North", 'short' = "South"))+
  # scale_y_continuous(limits = c(0,1.3)) +
  # geom_line(data = survreg_means_df, aes(color = core_edge, group = population), position = position_jitter(width = 0.2, seed = 5), show.legend = F) +
  geom_pointrange(data = survreg_means_df,
                  aes(x = treatment, y = emmean, ymax = emmean + SE, ymin = emmean - SE, color = core_edge),
                  position = position_jitter(width = 0.2, seed = 5), size = 0.9, alpha = 0.5, show.legend = T) +
  theme_classic2() +
  theme(text = element_text(size = 25), axis.ticks.x = element_blank())

# # Plot for local adaptation with lines ##Old way - with dodge
# ggplot(data = survreg_means_df, aes(x = treatment_core_edge, y = emmean)) +
#   geom_pointrange(aes(ymax = emmean + SE, ymin = emmean - SE, color = core_edge, group = population),
#                   position = position_dodge(width = 0.5), size = 0.9, alpha = 0.5, lty = 'solid') +
#   geom_line(aes(color = core_edge, group = population), position = position_dodge(width = 0.5), show.legend = F) +
#   geom_pointrange(data = survreg_means_CE_df, aes(x = treatment_core_edge, y = emmean, ymax = emmean + SE, ymin = emmean - SE, color = core_edge),
#              position = position_dodge(width = -0.5), stroke = 2, size = 1.3, show.legend = F) +
#   scale_color_manual(values = c(core_col, edge_col, NA), labels = c('core' = 'North', 'edge' = 'South')) +
#   scale_x_discrete(breaks = c('long_core', 'short_core'), labels = c("North", "South"))+
#   labs(y = "Days until Diapause", x = "Fall Environment", color = "Collection \nLocation") +
#   theme_classic() +
#   theme(text = element_text(size = 25), axis.ticks.x = element_blank(),
#         axis.text.x=element_text(hjust=-0.4))
# ggsave(file="lines.svg", path = 'G://My Drive//GRAD SCHOOL//CONFERENCES//Evolution 2022', width=7.53, height=5.4, units = 'in')

# # Plot for local adaptation without lines - doesn't work
# ggplot(data = survreg_means_df, aes(x = treatment, y = emmean, group = core_edge)) +
#   geom_errorbar(aes(ymax = emmean + SE, ymin = emmean - SE), position = position_dodge(width = 0.5), width = 0) +
#   geom_point(aes(color = core_edge), position = position_dodge(width = 0.5), size = 4, alpha = 0.5) +
#   geom_line(aes(color = core_edge, group = interaction(core_edge, population)), position = position_dodge2(width = 0.5)) +
#   geom_errorbar(data = survreg_means_CE_df,
#                 aes(ymax = emmean + SE, ymin = emmean - SE), position = position_dodge(width = -0.5), width = 0.3) +
#   geom_point(data = survreg_means_CE_df, aes(x = treatment, y = emmean, color = core_edge),
#              position = position_dodge(width = -0.5), pch = 1, stroke = 2, size = 5) +
#   scale_color_manual(values = c(core_col, edge_col), labels = c('core' = 'North', 'edge' = 'South')) +
#   labs(y = "Days until Diapause", x = "Fall Environment", color = "Collection \nLocation") +
#   scale_x_discrete(labels = c("long" = "North", "short" = "South")) +
#   theme_classic() +
#   theme(text = element_text(size = 17))

# Survival Regression Models with Latitude---- 
## Population
model1c_2 <- survreg(survObj ~ core_edge * latitude * treatment, dist = "gaussian", data = diapauseData3)
summary(model1c_2)
anova(model1c_2)
emtrends(model1c_2, specs = ~  treatment | core_edge, var = "latitude") ## results here
emmeans(model1c_2, pairwise ~ core_edge | treatment, type = "response")

# Other survival models ----
model1b <- survfit(survObj ~ core_edge + treatment, data = diapauseData3)
model1b
summary(model1b)

ggsurvplot(model1b, data = diapauseData3, 
           conf.int = T,
           surv.median.line = "hv",
           facet.by = "treatment",
           legend.title = "population",
           legend.labs = c("north", "south"),
           ylab = "Probability not in diapause",
           xlab = "Time (days)",
           censor.shape = 124,
           censor = F,
           ggtheme = theme_bw(),
           legend = "right")

#cox ph model - core_edge
model2 <- coxph(survObj ~ treatment * core_edge, data = diapauseData3)
model2
summary(model2)
anova(model2)
emmeans(model2, pairwise ~ core_edge|treatment, type = "surv") #type = "response")

#cox ph model - populations
model2b <- coxph(survObj ~ treatment + population, data = diapauseData3)
model2b
summary(model2b)
anova(model2b)
emmeans(model2b, pairwise ~ treatment|population, type = "response") #, back.transform = T) #  #type = "surv"

model2b_means <- emmeans(model2b, pairwise ~ treatment|population, type = "response")
model2b_means_df <- as.data.frame(model2b_means$emmeans) %>%
  mutate(core_edge = case_when(
    population == "Bl" | population == "Lc" | population == "Lj" | population == "Wi" ~ "edge",
    population == "De" | population == "Hu" | population == "Lo" | population == "Pu" ~ "core"
  ))

ggplot(data = model2b_means_df, aes(x = population, y = hazard)) +
  geom_errorbar(aes(ymax = asymp.UCL, ymin = asymp.LCL)) +
  geom_point(aes(color = treatment), size = 4) +
  facet_grid(~core_edge) +
  theme_bw()

# REGRESSION MODELS ----
## Diapausers only----
## dataset
diapauseDataDiapausing <- diapauseData3 %>%
  filter(diapause == 1)

## summarize data
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

##model with core_Edge, pop as rand
modDiapausing <- glmmTMB(time ~ core_edge * treatment + (1| core_edge / population), data = diapauseDataDiapausing,
                         family = poisson)
simulateResiduals(modDiapausing) %>% plot()
summary(modDiapausing)
Anova(modDiapausing, type = 3)
emmeans(modDiapausing, pairwise ~ core_edge|treatment, type = 'response')
emmeans(modDiapausing, pairwise ~ treatment|core_edge, type = 'response')

ranef(modDiapausing, condVar = F)
coef(modDiapausing)
ggpredict(modDiapausing, terms = c("population"), type = "random", condition = c(treatment = 'long'))

emoutDiapausing <- emmeans(modDiapausing, pairwise ~ core_edge*treatment, type = "response")
meansDiapausing_df <- as.data.frame(emoutDiapausing$emmeans)
meansDiapausing_df$core_edge <- factor(meansDiapausing_df$core_edge, levels = c('core', 'edge'))
meansDiapausing_df <- cbind(meansDiapausing_df,dtd_n1)

diapauseDataDiapausing$core_edge <- factor(diapauseDataDiapausing$core_edge, levels = c('core', 'edge'))


# Reaction Norm Plot
diapausers_coreedge_plot <- ggplot(data = meansDiapausing_df, aes(x = treatment, y = rate, fill = core_edge)) +
  geom_point(data = diapauseDataDiapausing, aes(x = treatment, y = time, fill = core_edge),
             position = position_jitterdodge(dodge.width = 0.35, jitter.width = 0.15),
             alpha = 0.75, shape = 21) +
  # geom_violin(data = diapauseDataDiapausing, aes(x = treatment, y = time, fill = core_edge), 
  #             position = position_jitterdodge(dodge.width = 0.35, jitter.width = 0.05),
  #             alpha = 0.25) +
  geom_line(aes(group = core_edge), position = position_dodge(width = 0.35), size = 0.73) +
  geom_pointrange(aes(x = treatment, y = rate, ymin = lower.CL, ymax = upper.CL), 
                  position = position_dodge(width = 0.35), shape = 21, size = 0.8) +
  geom_text(aes(label = dtd_n1, y = -1), position = position_dodge(width = 0.5)) +
  labs(x = "Test Environment", y = "Days to Diapause (95% CI)", fill = "Collection \nLocation") +
  scale_x_discrete(labels = c('long' = "North Fall", 'short' = "South Fall")) +
  scale_fill_manual(values = c('core' = 'white', 'edge' = 'black'), labels = c('edge' = "South", 'core' = "North")) +
  theme_bw(base_size = 15)
diapausers_coreedge_plot

# Bar Plot
ggplot(data = meansDiapausing_df) +
  geom_col(aes(x = core_edge, y = rate), color = 'black', fill = 'white', size = 1) +
  geom_errorbar(aes(x = core_edge, ymin = lower.CL, ymax = upper.CL), width = 0.4, size = 1) +
  facet_grid( ~ treatment, labeller = labeller(treatment = c('long' = "North Fall Environment", 'short' = "South Fall Environment"))) +
  geom_point(data = diapauseDataDiapausing, aes(x = core_edge, y = time), 
             position = position_jitter(width = 0.25), alpha = 0.25) +
  labs(x = "Collection Location", y = "Days to Diapause (95% CI)", fill = "Collection \nLocation") +
  scale_x_discrete(labels = c('edge' = "South", 'core' = "North")) +
  theme_bw(base_size = 15) +
  theme(strip.background = element_rect(fill = 'white'))

# Bar Plot #2
ggplot(data = meansDiapausing_df) +
  geom_col(aes(x = treatment, y = rate), color = 'black', fill = 'white', size = 1) +
  geom_errorbar(aes(x = treatment, ymin = lower.CL, ymax = upper.CL), width = 0.4, size = 1) +
  facet_grid( ~ core_edge, labeller = labeller(core_edge = c('core' = "North Collections", 'edge' = "South Collections"))) +
  geom_point(data = diapauseDataDiapausing, aes(x = treatment, y = time), 
             position = position_jitter(width = 0.25), alpha = 0.25) +
  labs(x = "Test Environment", y = "Days to Diapause (95% CI)", fill = "Collection \nLocation") +
  scale_x_discrete(labels = c('long' = "North Fall", 'short' = "South Fall ")) +
  theme_bw(base_size = 15) +
  theme(strip.background = element_rect(fill = 'white'))

#Plot of population raw means
diapausingPopMeansRaw <- diapauseDataDiapausing %>% group_by(treatment, population) %>% summarise(
  mean = mean(time),
  sd = sd(time),
  n = n(),
  se = sd/sqrt(n)
)
ggplot(data = diapausingPopMeansRaw) +
  geom_col(aes(x = population, y = mean), color = 'black', fill = 'white', size = 1) +
  geom_errorbar(aes(x = population, ymin = mean-se, ymax = mean+se), width = 0.4, size = 1) +
  facet_grid( ~ treatment, labeller = labeller(treatment = c('long' = "North Fall Environment", 'short' = "South Fall Environment"))) +
  geom_point(data = diapauseDataDiapausing, aes(x = population, y = time),
             position = position_jitter(width = 0.25), alpha = 0.25) +
  labs(title = 'Population Raw Means and Standard Errors', x = "Collection Location", y = "Days to Diapause (SE)", fill = "Collection \nLocation") +
  # scale_x_discrete(labels = c('edge' = "South", 'core' = "North")) +
  theme_bw(base_size = 15) +
  theme(strip.background = element_rect(fill = 'white'))
  
##model with populations - low sample sizes here!
modDiapausingpop <- glmmTMB(time ~  treatment * population, data = diapauseDataDiapausing,
                            family = poisson)
simulateResiduals(modDiapausingpop) %>% plot()
summary(modDiapausingpop)
Anova(modDiapausingpop, type = 3)
emmeans(modDiapausingpop, pairwise ~ treatment*population, type = 'response')

emoutDiapausingPops <- emmeans(modDiapausingpop, pairwise ~ population*treatment, type = "response")
meansDiapausingPops_df <- as.data.frame(emoutDiapausingPops$emmeans) %>%
  mutate(core_edge = case_when(
    population == "Bl" | population == "Lc" | population == "Lj" | population == "Wi" ~ "edge",
    population == "De" | population == "Hu" | population == "Lo" | population == "Pu" ~ "core"
  ))
meansDiapausingPops_df$population <- factor(meansDiapausingPops_df$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))
diapauseDataDiapausing$population <- factor(diapauseDataDiapausing$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))

diapausers_pop_plot <- ggplot(data = meansDiapausingPops_df) +
  geom_point(data = diapauseDataDiapausing, aes(x = population, y = time),
             position = position_jitter(width = 0.15), alpha = 0.25) +
  geom_pointrange(aes(x = population, y = rate, ymin = lower.CL, ymax = upper.CL, fill = core_edge), 
                  color = 'black', shape = 21, size = 0.8) +
  facet_grid( ~ treatment, labeller = labeller(treatment = c('long' = "North Fall Environment", 'short' = "South Fall Environment")),
              switch = 'x') +
  geom_text(data = dtd_n2, aes(label = n2, x = population, y = -1)) +
  labs(x = "Collection Location (N to S)", y = "Days to Diapause (95% CI)", fill = "Collection \nLocation") +
  scale_x_discrete(labels = c('Lo' = 'A', 'Hu' = 'B', 'De' = 'C', 'Pu' = 'D', 'Lc' = 'E', 'Wi' = 'F', 'Lj' = 'G', 'Bl' = 'H')) +
  scale_fill_manual(values = c('core' = 'white', 'edge' = 'black'), labels = c('edge' = "South", 'core' = "North")) +
  theme_bw(base_size = 15) +
  theme(strip.background = element_rect(fill = 'white'))
diapausers_pop_plot


##Combining plots for publication
ggarrange(diapausers_coreedge_plot, diapausers_pop_plot, 
          common.legend = T, 
          legend = 'right', 
          widths = c(1.3,2),
          labels = 'AUTO')
##export: width = 1100, height = 575 (11.46in x 5.99in)


##model with populations for only core
modDiapausingpopcore <- glmmTMB(time ~  treatment * population, data = diapauseDataDiapausing %>% filter(core_edge == 'core'),
                            family = poisson)
simulateResiduals(modDiapausingpopcore) %>% plot()
summary(modDiapausingpopcore)
Anova(modDiapausingpopcore, type = 3)
emmeans(modDiapausingpopcore, pairwise ~ treatment|population, type = "response")

emoutDiapausingCore <- emmeans(modDiapausingpopcore, pairwise ~ population*treatment, type = "response")
meansDiapausingCore_df <- as.data.frame(emoutDiapausingCore$emmeans)
# meansDiapausingCore_df$population <- factor(meansDiapausingCore_df$population, levels = c('Lo', 'Hu', 'De', 'Pu'))


##model with populations for only edge
modDiapausingpopedge <- glmmTMB(time ~  treatment * population, data = diapauseDataDiapausing %>% filter(core_edge == 'edge'),
                                family = poisson)
simulateResiduals(modDiapausingpopedge) %>% plot()
summary(modDiapausingpopedge)
Anova(modDiapausingpopedge, type = 3)
emmeans(modDiapausingpopedge, pairwise ~ treatment|population, type = "response")

emoutDiapausingEdge <- emmeans(modDiapausingpopedge, pairwise ~ population*treatment, type = "response")
meansDiapausingEdge_df <- as.data.frame(emoutDiapausingEdge$emmeans)
# meansDiapausingEdge_df$population <- factor(meansDiapausingEdge_df$population, levels = c('Lc', 'Wi', 'Lj', 'Bl'))

##Plot for populations, from separate models
diapauseDataDiapausing$population <- factor(diapauseDataDiapausing$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))

meansDiapausingPopsBind_df <- rbind(meansDiapausingCore_df, meansDiapausingEdge_df)
meansDiapausingPopsBind_df$population <- factor(meansDiapausingPopsBind_df$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))


ggplot(data = meansDiapausingPopsBind_df) +
  geom_col(aes(x = population, y = rate), color = 'black', fill = 'white', size = 1) +
  geom_errorbar(aes(x = population, ymin = lower.CL, ymax = upper.CL), width = 0.4, size = 1) +
  facet_grid( ~ treatment, labeller = labeller(treatment = c('long' = "North Fall Environment", 'short' = "South Fall Environment"))) +
  geom_point(data = diapauseDataDiapausing, aes(x = population, y = time),
             position = position_jitter(width = 0.25), alpha = 0.25) +
  labs(title = 'Population Means (Separate Models) and 95% CI', x = "Collection Location", y = "Days to Diapause (95% CI)", fill = "Collection \nLocation") +
  # scale_x_discrete(labels = c('edge' = "South", 'core' = "North")) +
  theme_bw(base_size = 15) +
  theme(strip.background = element_rect(fill = 'white'))

#Models with Latitude
modDiapausinglat <- glmmTMB(time ~  treatment * core_edge * latitude, data = diapauseDataDiapausing,
                                family = poisson)
simulateResiduals(modDiapausinglat) %>% plot()
summary(modDiapausinglat)
Anova(modDiapausinglat, type = 3)
emtrends(modDiapausinglat, specs = ~treatment*core_edge, var = "latitude", regrid = 'response')
emtrends(modDiapausinglat, specs = ~treatment*core_edge, var = "latitude", type = 'response')

emmip(modDiapausinglat, treatment | core_edge ~ latitude, cov.reduce = range, CIs = F)
emmeans(modDiapausinglat, pairwise ~ treatment*core_edge*latitude, type = 'response', at = list(latitude = c(43, 38, 33)))

#at = list(latitude = c(44.856, 40.063, 39.144, 38.268, 34.593, 34.422, 34.342, 33.912))


# ##model with latitude Core
# modDiapausinglatCore <- glmmTMB(time ~  treatment * latitude, data = diapauseDataDiapausing %>% filter(core_edge == 'core') ,
#                                 family = poisson)
# simulateResiduals(modDiapausinglatCore) %>% plot()
# summary(modDiapausinglatCore)
# Anova(modDiapausinglatCore, type = 3)
# emtrends(modDiapausinglatCore, specs = ~treatment, var = "latitude", type = 'response')
# 
# ##model with latitude Edge
# 
# modDiapausinglatEdge <- glmmTMB(time ~  treatment * latitude, data = diapauseDataDiapausing %>% filter(core_edge == 'edge') ,
#                             family = poisson)
# simulateResiduals(modDiapausinglatEdge) %>% plot()
# summary(modDiapausinglatEdge)
# Anova(modDiapausinglatEdge, type = 3)
# emtrends(modDiapausinglatEdge, specs = ~treatment, var = "latitude", type = 'response')


  

## All samples ----

# create new column where non diapausers are 43, diapausers are days to diapause
diapauseData3 <- diapauseData3 %>% mutate(time_diap = replace(time, diapause == 0, 43))
diapauseData3$core_edge <- factor(diapauseData3$core_edge, levels = c('core', 'edge'))
diapauseData3$degree_days_scaled <- scale(diapauseData3$degree_days, center = T, scale = T)
diapauseData3$population <- factor(diapauseData3$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))


## summarize data
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


ggplot(data = diapauseData3, aes(color = interaction(treatment, core_edge))) +
  geom_density(aes(x = time))

library(multimode)
modetest(diapauseData3$time)
locmodes(diapauseData3$time, mod0 = 2, display = TRUE)

## model with core_Edge, pop as rand
modAllSamp <- glmmTMB(time_diap ~ core_edge * treatment + (1| population), data = diapauseData3,
                         family = poisson())
simulateResiduals(modAllSamp) %>% plot()
summary(modAllSamp)
Anova(modAllSamp, type = 3)
emmeans(modAllSamp, pairwise ~ core_edge|treatment, type = 'response')
emmeans(modAllSamp, pairwise ~ treatment|core_edge, type = 'response')

ranef(modAllSamp, condVar = F)
coef(modAllSamp)
ggpredict(modAllSamp, terms = c("population"), type = "random", condition = c(treatment = 'long'))

emoutmodAllSamp <- emmeans(modAllSamp, pairwise ~ core_edge*treatment, type = "response")
emoutmodAllSamp_df <- as.data.frame(emoutmodAllSamp$emmeans)
emoutmodAllSamp_df$core_edge <- factor(emoutmodAllSamp_df$core_edge, levels = c('core', 'edge'))


# Reaction Norm Plot
allSamp_coreedge_plot <- ggplot(data = emoutmodAllSamp_df, aes(x = treatment, y = rate, fill = core_edge)) +
  annotate('rect', xmin = c(0.8,1.8), xmax = c(1.2, 2.2), ymin = 42, ymax = 44, alpha = 0.2, color = "black") +
  geom_point(data = diapauseData3, aes(x = treatment, y = time_diap, fill = core_edge),
             position = position_jitterdodge(dodge.width = 0.35, jitter.width = 0.15),
             alpha = 0.25) +
  geom_line(aes(group = core_edge), position = position_dodge(width = 0.35), size = 0.73) +
  geom_pointrange(aes(x = treatment, y = rate, ymin = lower.CL, ymax = upper.CL),
                  position = position_dodge(width = 0.35), shape = 21, size = 0.8) +
  geom_text(data = summ_dtd_all, aes(label = n2, y = -1.5), position = position_dodge(width = 0.35)) +
  annotate("label", x = c(1,2), y = 45.5, label = "Non-Diapausers") +
  labs(x = "Diapause-Inducing Daylength", y = "Days to Diapause (95% CI)", fill = "Collection \nLocation") +
  scale_x_discrete(labels = c('long' = "North", 'short' = "South")) +
  scale_fill_manual(values = c('core' = 'white', 'edge' = 'black'), labels = c('edge' = "South", 'core' = "North")) +
  theme_bw(base_size = 15)
allSamp_coreedge_plot



# Population analyses - all samples
modAllSamppop <- glmmTMB(time_diap ~  treatment * population, data = diapauseData3,
                            family = poisson())
simulateResiduals(modAllSamppop) %>% plot()
summary(modAllSamppop)
Anova(modAllSamppop, type = 3)
emmeans(modAllSamppop, pairwise ~ treatment|population, type = 'response')

emoutAllSampPops <- emmeans(modAllSamppop, pairwise ~ population*treatment, type = "response")
emoutAllSampPops_df <- as.data.frame(emoutAllSampPops$emmeans) %>%
  mutate(core_edge = case_when(
    population == "Bl" | population == "Lc" | population == "Lj" | population == "Wi" ~ "edge",
    population == "De" | population == "Hu" | population == "Lo" | population == "Pu" ~ "core"
  ))
emoutAllSampPops_df$population <- factor(emoutAllSampPops_df$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))

allSamp_pop_plot <- ggplot(data = emoutAllSampPops_df) +
  annotate('rect', xmin = 0.75, xmax = 8.25, ymin = 42, ymax = 44, alpha = 0.2, color = 'black') +
  geom_point(data = diapauseData3, aes(x = population, y = time_diap),
             position = position_jitter(width = 0.15), alpha = 0.25) +
  geom_pointrange(aes(x = population, y = rate, ymin = lower.CL, ymax = upper.CL, fill = core_edge),
                  color = 'black', shape = 21, size = 0.8) +
  facet_grid( ~ treatment, labeller = labeller(treatment = c('long' = "Northern Daylength", 'short' = "Southern Daylength"))) +
  geom_text(data = summ_dtd_all_pop, aes(label = n2, x = population, y = -1)) +
  annotate("label", x = 4.5, y = 46, label = "Non-Diapausers") +
  labs(x = "Collection Location (N to S)", y = "Days to Diapause (95% CI)", fill = "Collection \nLocation") +
  scale_x_discrete(labels = c('Lo' = 'A', 'Hu' = 'B', 'De' = 'C', 'Pu' = 'D', 'Lc' = 'E', 'Wi' = 'F', 'Lj' = 'G', 'Bl' = 'H')) +
  scale_fill_manual(values = c('core' = 'white', 'edge' = 'black'), labels = c('edge' = "South", 'core' = "North")) +
  theme_bw(base_size = 15) +
  theme(strip.background = element_rect(fill = 'white'))
allSamp_pop_plot

##Combining figures for publication
ggarrange(allSamp_coreedge_plot, allSamp_pop_plot, 
          common.legend = T, 
          legend = 'right', 
          widths = c(1.3,2),
          labels = 'AUTO')
##export: width = 1100, height = 575 (11.46in x 5.99in)

## models with latitude - all samples
modAllSamplat <- glmmTMB(time_diap ~  treatment * core_edge * latitude, data = diapauseData3,
                            family = poisson(link = 'log'))
simulateResiduals(modAllSamplat) %>% plot()
summary(modAllSamplat)
Anova(modAllSamplat, type = 3)
emtrends(modAllSamplat, specs = ~treatment*core_edge, var = "latitude", regrid = "response")

emmip(modAllSamplat, treatment | core_edge ~ latitude, cov.reduce = range, regrid = "response", CIs = T)
emmeans(modAllSamplat, pairwise ~ treatment*core_edge*latitude, type = 'response', at = list(latitude = c(43, 38, 33)))


## models with other site characteristics - all samples
modAllSamp_sitechar <- glmmTMB(time_diap ~  (treatment + core_edge + scale(latitude) + degree_days_scaled)^3, data = diapauseData3,
                         family = gaussian)
                               # family = poisson(link = 'log'))
simulateResiduals(modAllSamp_sitechar) %>% plot()
summary(modAllSamp_sitechar)
Anova(modAllSamp_sitechar, type = 3)
emtrends(modAllSamp_sitechar, specs = ~treatment*core_edge, var = "latitude", regrid = "response")


modAllSamp_edgeDD <- glmmTMB(time_diap ~  latitude, 
                             data = diapauseData3 %>% filter(core_edge == 'edge') %>% filter(treatment == 'short'),
                               family = gaussian)
simulateResiduals(modAllSamp_edgeDD) %>% plot()
summary(modAllSamp_edgeDD)
Anova(modAllSamp_edgeDD, type = 3)
emtrends(modAllSamp_edgeDD, var = "elevation", regrid = "response")

cor.test(diapauseData3$time_diap, diapauseData3$degree_days)
cor.test(site_char$elevation, site_char$degree_days)

##which of the three variables best predicts responses in core and edge separately - based on r2 or AIC??
## Edge, Short 
e_s <- diapauseData3 %>% filter(core_edge == 'edge') %>% filter(treatment == 'short')

mod_e_s_lat <- lm(time_diap ~ latitude, data = e_s)
summary(mod_e_s_lat)

mod_e_s_ele <- lm(time_diap ~ elevation, data = e_s)
summary(mod_e_s_ele)

mod_e_s_deg <- lm(time_diap ~ degree_days, data = e_s)
summary(mod_e_s_deg)

#Plots

site_char$letter <- c("Site H", "Site E", "Site G", "Site F", "Site C", "Site B", "Site A", "Site D")

e_lat <- ggplot(data = e_s, aes(x = latitude, y = time_diap)) +
  geom_smooth(method = lm) +
  geom_jitter(width = .01, height = 2) +
  theme_bw(base_size = 15) +
  labs(y = "Days to Diapause", x = "Latitude") +
  annotate("text", y = 42, x = 34.1, label = "italic(R) ^ 2 == 0.226", parse = TRUE, size = 8) +
  annotate("text", y = -2, x = site_char$latitude[1:4], label = site_char$letter[1:4], size = 5)
e_lat

e_ele <- ggplot(data = e_s, aes(x = elevation, y = time_diap)) +
  geom_smooth(method = lm) +
  geom_jitter(width = 15, height = 2) +
  theme_bw(base_size = 15) +
  labs(y = "Days to Diapause", x = "Elevation") +
  annotate("text", y = 42, x = 500, label = "italic(R) ^ 2 == 0.340", parse = TRUE, size = 8) +
  annotate("text", y = -2, x = site_char$elevation[1:4], label = site_char$letter[1:4], size = 5)
e_ele

e_deg <- ggplot(data = e_s, aes(x = degree_days, y = time_diap)) +
  geom_smooth(method = lm) +
  geom_jitter(width = 15, height = 2) +
  theme_bw(base_size = 15) +
  labs(y = "Days to Diapause", x = "Cumulative Degree Days") +
  annotate("text", y = 42, x = 2500, label = "italic(R) ^ 2 == 0.494", parse = TRUE, size = 8) +
  annotate("text", y = -2, x = site_char$degree_days[1:4], label = site_char$letter[1:4], size = 5)
e_deg

ggarrange(e_lat,e_ele,e_deg, nrow = 1)
#export: 1750 x 700 px

## Core, Long 
c_l <- diapauseData3 %>% filter(core_edge == 'core') %>% filter(treatment == 'long')

mod_c_l_lat <- lm(time_diap ~ latitude, data = c_l)
summary(mod_c_l_lat)

mod_c_l_ele <- lm(time_diap ~ elevation, data = c_l)
summary(mod_c_l_ele)

mod_c_l_deg <- lm(time_diap ~ degree_days, data = c_l)
summary(mod_c_l_deg)

##
mod_c_l_all <- lm(time_diap ~ latitude+elevation+degree_days, data = c_l)
summary(mod_c_l_all)
Anova(mod_c_l_all, type = 3)

# modAllSamplat_glm <- glm(time_diap ~  treatment * core_edge * latitude, data = diapauseData3,
#      family = poisson(link = 'log'))
# summary(modAllSamplat_glm)
# Anova(modAllSamplat_glm, type = 3)
# emtrends(modAllSamplat_glm, specs = ~treatment*core_edge, var = "latitude", regrid = "response")




# PROPORTION DIAPAUSE ANALYSES ----

## Logistic Regression

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

diapauseData3$population <- factor(diapauseData3$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))
diapauseData3$core_edge <- factor(diapauseData3$core_edge, levels = c('core', 'edge'))

diapauseData3Aug <- diapauseData3 %>%   #add_row(treatment = "short", core_edge = "core", diapause = 0)
  add_row(treatment = "short", core_edge = "core", population = 'Hu', diapause = 0) %>%
  add_row(treatment = "short", core_edge = "core", population = 'Lo', diapause = 0) %>%
  add_row(treatment = "short", core_edge = "core", population = 'Pu', diapause = 0) %>%
  add_row(treatment = "short", core_edge = "edge", population = 'Lj', diapause = 0)
diapauseData3Aug$population <- factor(diapauseData3Aug$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))
diapauseData3Aug$core_edge <- factor(diapauseData3Aug$core_edge)
diapauseData3Aug$treatment <- factor(diapauseData3Aug$treatment)

diapauseData3Trim <- diapauseData3 %>% 
  filter(population != 'Hu' | treatment != 'short') %>% 
  filter(population != 'Lo' | treatment != 'short') %>% 
  filter(population != 'Pu' | treatment != 'short') %>% 
  filter(population != 'Lj' | treatment != 'short')
  

logmodel1 <- glmmTMB(diapause ~ core_edge * treatment + (1|population), family = "binomial", data = diapauseData3)
simulateResiduals(logmodel1) %>% plot()
summary(logmodel1)
Anova(logmodel1, type = 3)
ranef(logmodel1, condVar = F)
coef(logmodel1)
emmeans(logmodel1, pairwise ~ core_edge | treatment, type = "response")

logistic_means_ce <- emmeans(logmodel1, pairwise ~ treatment * core_edge, type = "response")
logistic_means_ce_df <- as.data.frame(logistic_means_ce$emmeans) %>%
  unite("treatment_core_edge", c(treatment,core_edge), remove = F)

# Reaction Norm Plot
prop_coreedge_plot <- ggplot(data = logistic_means_ce_df, aes(x = treatment, y = prob, fill = core_edge)) +
  geom_line(aes(group = core_edge), position = position_dodge(width = 0.35), size = 0.73) +
  # geom_point(data = diapauseData3, aes(x = treatment, y = diapause, fill = core_edge),
  #            position = position_jitterdodge(dodge.width = 0.35, jitter.width = 0.15, jitter.height = 0.05),
  #            alpha = 0.25, shape = 16) +
  geom_pointrange(aes(x = treatment, y = prob, ymin = lower.CL, ymax = upper.CL), 
                  position = position_dodge(width = 0.35), shape = 21, size = 0.8, color = 'black') +
  labs(x = "Diapause-Inducing Daylength", y = "Proportion in Diapause (95% CI)", fill = "Collection \nLocation") +
  scale_x_discrete(labels = c('long' = "North", 'short' = "South")) +
  scale_fill_manual(values = c('core' = 'white', 'edge' = 'black'), labels = c('edge' = "South", 'core' = "North")) +
  theme_bw(base_size = 15)
prop_coreedge_plot

# Bar Plot
ggplot(data = logistic_means_ce_df) +
  geom_col(aes(x = core_edge, y = prob), color = 'black', fill = 'white', size = 1) +
  geom_errorbar(aes(x = core_edge, ymin = lower.CL, ymax = upper.CL), width = 0.4, size = 1) +
  facet_grid( ~ treatment, labeller = labeller(treatment = c('long' = "North Fall Environment", 'short' = "South Fall Environment"))) +
  geom_point(data = diapauseData3, aes(x = core_edge, y = diapause), 
             position = position_jitter(height = 0.05, width = 0.35), alpha = 0.25) +
  labs(x = "Collection Location", y = "Proportion in Diapause (95% CI)", fill = "Collection \nLocation") +
  scale_x_discrete(labels = c('edge' = "South", 'core' = "North")) +
  theme_bw(base_size = 15) +
  theme(strip.background = element_rect(fill = 'white'))

##latitude model
logmodel2 <- glmmTMB(diapause ~ latitude * treatment * core_edge, family = binomial , data = diapauseData3Aug)
simulateResiduals(logmodel2) %>% plot()
summary(logmodel2)
Anova(logmodel2, type = 3)
emtrends(logmodel2, specs = ~treatment*core_edge, var = "latitude", regrid = 'response')
# emtrends(logmodel2, specs = ~treatment*core_edge, var = "latitude", type = 'response')
# emtrends(logmodel2, specs = ~treatment*core_edge, var = "latitude")
# 
# emtrends(logmodel2, pairwise ~treatment*core_edge, var = "latitude", regrid = 'response')

emmip(logmodel2, treatment | core_edge ~ latitude, cov.reduce = range, regrid = "response", CI = T)


ggplot(data = diapusePropSummaryPop, aes(x = latitude, y = prop, color = population)) +
  geom_point(aes(shape = treatment)) +
  facet_grid(~treatment) +
  # geom_line(aes(group = population), position = position_dodge(width = 0.35), size = 0.73) +
  # geom_point(data = diapauseData3, aes(x = treatment, y = diapause, fill = core_edge),
  #            position = position_jitterdodge(dodge.width = 0.35, jitter.width = 0.15, jitter.height = 0.05),
  #            alpha = 0.25, shape = 16) +
  # geom_pointrange(aes(x = treatment, y = prob, ymin = lower.CL, ymax = upper.CL), 
  #                 position = position_dodge(width = 0.35), shape = 21, size = 0.8, color = 'black') +
  # labs(x = "Test Environment", y = "Proportion in Diapause (95% CI)", fill = "Collection \nLocation") +
  # scale_x_discrete(labels = c('long' = "North Fall", 'short' = "South Fall")) +
  # scale_fill_manual(values = c('core' = 'white', 'edge' = 'black'), labels = c('edge' = "South", 'core' = "North")) +
  theme_bw(base_size = 15)

#Firth penalized regression for population
library(logistf)
logmodelFirth <- logistf(diapause ~ population * treatment, pl = TRUE, data = diapauseData3)
summary(logmodelFirth)
logmodelFirth2 <- logistf(diapause ~ population + treatment, pl = TRUE, data = diapauseData3)
logmodelFirth3 <- logistf(diapause ~ treatment, pl = TRUE, data = diapauseData3)
logmodelFirth4 <- logistf(diapause ~ population, pl = TRUE, data = diapauseData3)

anova(logmodelFirth, logmodelFirth4)
tidy(coef(logmodelFirth))
emmeans(logmodelFirth, pairwise ~ treatment | population, type = "response")


#Firth regression with latitude
logmodel_pop_Firth <- logistf(diapause ~ latitude * treatment * core_edge, pl = TRUE, data = diapauseData3)



#
logmodel3Trim <- glmmTMB(diapause ~ population * treatment, family = 'binomial', data = diapauseData3Trim)


# logistic regression for population - augmented dataset
logmodel3 <- glmmTMB(diapause ~ population * treatment, family = 'binomial', data = diapauseData3Aug)
simulateResiduals(logmodel3) %>% plot()
summary(logmodel3)
Anova(logmodel3, type = 3)
emmeans(logmodel3, pairwise ~ treatment | population, type = "response")
logistic_means_pop <- emmeans(logmodel3, pairwise ~ treatment * population, type = "response")
logistic_means_pop_df <- as.data.frame(logistic_means_pop$emmeans) %>%
  mutate(core_edge = case_when(
    population == "Bl" | population == "Lc" | population == "Lj" | population == "Wi" ~ "edge",
    population == "De" | population == "Hu" | population == "Lo" | population == "Pu" ~ "core"
  ))
logistic_means_pop_df$population <- factor(logistic_means_pop_df$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))

# Plot of population proportions
prop_population_plot <- ggplot(data = logistic_means_pop_df) +
  geom_pointrange(aes(x = population, y = prob, ymin = lower.CL, ymax = upper.CL, fill = core_edge), 
                  shape = 21, size = 0.8) +
  # geom_col(aes(x = population, y = prob), color = 'black', fill = 'white', size = 1) +
  # geom_errorbar(aes(x = population, ymin = lower.CL, ymax = upper.CL), width = 0.4, size = 1) +
  facet_grid( ~ treatment, labeller = labeller(treatment = c('long' = "Northern Daylength", 'short' = "Southern Daylength"))) +
  # geom_point(data = diapauseData3, aes(x = population, y = diapause),
  #            position = position_jitter(height = 0.05, width = 0.25), alpha = 0.25) +
  labs(x = "Collection Location (N to S)", y = "Proportion in Diapause (95% CI)", fill = "Collection \nLocation") +
  scale_fill_manual(values = c('core' = 'white', 'edge' = 'black'), labels = c('edge' = "South", 'core' = "North")) +
    scale_x_discrete(labels = c('Lo' = 'A', 'Hu' = 'B', 'De' = 'C', 'Pu' = 'D', 'Lc' = 'E', 'Wi' = 'F', 'Lj' = 'G', 'Bl' = 'H')) +
  theme_bw(base_size = 15) +
  # ylim(-0.2,1.25) +
  theme(strip.background = element_rect(fill = 'white'))
prop_population_plot

##Combining plots for publication
ggarrange(prop_coreedge_plot, prop_population_bayes_plot, 
          common.legend = T, 
          legend = 'right', 
          widths = c(1.3,2),
          labels = 'AUTO')
##export: width = 1100, height = 575 (11.46in x 5.99in)

##Sample size 
prop_sample_size <- diapauseData3 %>% group_by(population, treatment) %>% summarize (
  n = n()
)
mean(prop_sample_size$n)

## Bayesian approach to fix complete separation ----
library(arm)
logmodel_bayes <- bayesglm(diapause ~ population * treatment, family = 'binomial', data = diapauseData3)
summary(logmodel_bayes)
Anova(logmodel_bayes, type = 3)
emmeans(logmodel_bayes, pairwise ~ treatment | population, type = "response")
logisticbayes_means_pop <- emmeans(logmodel_bayes, pairwise ~ treatment * population, type = "response")
logisticbayes_means_pop_df <- as.data.frame(logisticbayes_means_pop$emmeans) %>%
  mutate(core_edge = case_when(
    population == "Bl" | population == "Lc" | population == "Lj" | population == "Wi" ~ "edge",
    population == "De" | population == "Hu" | population == "Lo" | population == "Pu" ~ "core"
  ))
logisticbayes_means_pop_df$population <- factor(logisticbayes_means_pop_df$population, levels = c('Lo', 'Hu', 'De', 'Pu', 'Lc', 'Wi', 'Lj', 'Bl'))

#plot of bayes logistic approach
prop_population_bayes_plot <- ggplot(data = logisticbayes_means_pop_df) +
  geom_pointrange(aes(x = population, y = prob, ymin = asymp.LCL, ymax = asymp.UCL, fill = core_edge), 
                  shape = 21, size = 0.8) +
  # geom_col(aes(x = population, y = prob), color = 'black', fill = 'white', size = 1) +
  # geom_errorbar(aes(x = population, ymin = lower.CL, ymax = upper.CL), width = 0.4, size = 1) +
  facet_grid( ~ treatment, labeller = labeller(treatment = c('long' = "Northern Daylength", 'short' = "Southern Daylength"))) +
  # geom_point(data = diapauseData3, aes(x = population, y = diapause),
  #            position = position_jitter(height = 0.05, width = 0.25), alpha = 0.25) +
  labs(x = "Collection Location (N to S)", y = "Proportion in Diapause (95% CI)", fill = "Collection \nLocation") +
  scale_fill_manual(values = c('core' = 'white', 'edge' = 'black'), labels = c('edge' = "South", 'core' = "North")) +
  scale_x_discrete(labels = c('Lo' = 'A', 'Hu' = 'B', 'De' = 'C', 'Pu' = 'D', 'Lc' = 'E', 'Wi' = 'F', 'Lj' = 'G', 'Bl' = 'H')) +
  theme_bw(base_size = 15) +
  # ylim(-0.2,1.25) +
  theme(strip.background = element_rect(fill = 'white'))
prop_population_bayes_plot

#compare augmented and bayes approaches
ggarrange(prop_population_plot, prop_population_bayes_plot, common.legend = T)

## Bayes approach with latitude
logmodel_lat_bayes <- bayesglm(diapause ~ latitude * treatment * core_edge, family = binomial , data = diapauseData3)
simulateResiduals(logmodel_lat_bayes) %>% plot()
summary(logmodel_lat_bayes)
Anova(logmodel_lat_bayes, type = 3)
emtrends(logmodel_lat_bayes, specs = ~treatment*core_edge, var = "latitude", regrid = 'response')
emmip(logmodel_lat_bayes, treatment | core_edge ~ latitude, cov.reduce = range, regrid = "response", CI = T)


#####

library(lme4)
library(lmerTest)
logmodel5 <- glmer(diapause ~ core_edge * treatment + weight + (1|population), family = "binomial", data = diapauseData3)
summary(logmodel5)


##Plot of proportion in diapause
ggplot(data = logistic_means_ce_df, aes(x = treatment, y = prob)) +
  geom_pointrange(aes(ymax = prob + SE, ymin = prob - SE, color = core_edge), 
                  size = 1.3, stroke = 2, lty = 'solid') +
  labs(y = "Proportion in Diapause", x = "Fall Environment", color = "Collection \nLocation") +
  scale_color_manual(values = c(edge_col, core_col), labels = c('core' = 'North', 'edge' = 'South')) +
  scale_x_discrete(labels = c('long' = "North", 'short' = "South"))+
  # scale_y_continuous(limits = c(0,1.3)) +
  geom_line(data = logistic_means_pop_df, aes(color = core_edge, group = population), position = position_jitter(width = 0.2, seed = 5), show.legend = F) +
  geom_pointrange(data = logistic_means_pop_df,
  aes(x = treatment, y = prob, ymax = prob + SE, ymin = prob - SE, color = core_edge),
  position = position_jitter(width = 0.2, seed = 5), shape = 17, size = 0.9, alpha = 0.5, show.legend = T) +
  theme_classic2() +
  theme(text = element_text(size = 25), axis.ticks.x = element_blank())


# Population means for modeling

#proportion in diapause
diapauseData3 %>% group_by(treatment, population) %>% summarise(
  mean = mean(diapause, na.rm = TRUE),
  var = var(diapause, na.rm = TRUE),
  n = sum(is.na(diapause) == F))

#days until diapause
diapauseData3 %>% group_by(treatment, population) %>% summarise(
  mean = mean(time, na.rm = TRUE),
  var = var(time, na.rm = TRUE),
  n = sum(is.na(time) == F))


##Distance from core analyses
logmodel4 <- glmmTMB(diapause ~ core_edge + distFromCore + weight + (1|population), family = "binomial", 
                     data = diapauseData3 %>% filter(treatment == "long"))
simulateResiduals(logmodel4) %>% plot()
summary(logmodel4)
Anova(logmodel4, type = 3)
emmeans(logmodel4, pairwise ~ core_edge, type = "response")
# logistic_means_ce <- emmeans(logmodel1, pairwise ~ treatment * core_edge, type = "response")
# logistic_means_ce_df <- as.data.frame(logistic_means_ce$emmeans) %>%
#   unite("treatment_core_edge", c(treatment,core_edge), remove = F)


#Trying to transform data ----
# Data where no diapause = 0 days

diapauseData3 <- diapauseData3 %>% 
  mutate(time_non1 = if_else(diapause == 0, 0, time)) %>%
  mutate(time_non2 = if_else(diapause == 0, 0, time + 42))
hist(diapauseData3$time_non1, breaks = 50)


daysmodelNon <- glmmTMB(time_non1 ~ core_edge * treatment + weight + (1|population), family = "poisson", 
                     zi = ~ 1, data = diapauseData3)
simulateResiduals(daysmodelNon) %>% plot()
summary(daysmodelNon)
Anova(daysmodelNon, type = 3)
emmeans(daysmodelNon, pairwise ~ core_edge | treatment, type = "response")

daysmodelNon2 <- glmmTMB(time_non1 ~ population * treatment + weight, family = "poisson", 
                        # zi = ~ 1, 
                        data = diapauseData3)
simulateResiduals(daysmodelNon2) %>% plot()
summary(daysmodelNon2)
Anova(daysmodelNon2, type = 3)
emmeans(daysmodelNon2, pairwise ~ treatment | population, type = "response")




