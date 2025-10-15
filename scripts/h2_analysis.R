##############################################################################.
############.
############         Estimate heritability and evolvability
############      of D. carinulata Days to Diapause
############.      
############                Written by: ############
############            Date Last Modified: 2/24/2025
############.      
##############################################################################.



# Packages ----------------------------------------------------------------

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggpubr)
library(merTools)

select <- dplyr::select


# Custom Functions ---------------------------------------------------------

## Extract variance components from a random effects model

var_from_model <- function(model, mean, design) {
  
  model_results <- as.data.frame(VarCorr(model))
  
  # For paternal half-sib design
  if(design == 'phs'){
    REML.results <- data.frame(var_sire = NA, var_e = NA, V_A = NA, V_P = NA, h2 = NA, I = NA) #blank data.frame to store results
    REML.results.all <- REML.results %>% add_row(var_sire = model_results[1,'vcov'], #store results and calculate VA & VP & h2
                                                 var_e = model_results[2,'vcov'],
                                                 V_A = 4 * var_sire,
                                                 V_P = var_sire + var_e,
                                                 h2 = V_A / V_P,
                                                 I = V_A / mean^2) %>%
      filter(is.na(var_sire) == F)
    
    REML.results2 <- REML.results.all %>%
      dplyr::select(c(V_A:I)) %>%
      pivot_longer(cols = c(1:4), names_to = "var_comp", values_to = "est")
    REML.results2$var_comp <- factor(REML.results2$var_comp, levels = c('V_P', 'V_A', 'h2', 'I')) #change order of VA and VP
    levels(REML.results2$var_comp) <- c('V[Phenotypic]', 'V[Additive]', 'h^2', 'I[A]') #change labels for VA and VP
    
    output = list(AllVarComps = REML.results.all, VarCompsLong = REML.results2)
    
    return(output)
  }
  
  # For paternal half-sib-full-sib design
  if(design == 'fs'){
    REML.results <- data.frame(var_sire = NA, var_dam = NA, var_e = NA, V_A = NA, V_P = NA, h2 = NA, I = NA) #blank data.frame to store results
    REML.results.all <- REML.results %>% add_row(var_sire = c(model_results[2,'vcov']), #store results and calculate VA & VP & h2
                                                 var_dam = c(model_results[1,'vcov']),
                                                 var_e = c(model_results[3,'vcov']), 
                                                 V_A = c(4 * var_sire),
                                                 V_P = c(var_sire + var_dam + var_e),
                                                 h2 = c(V_A / V_P),
                                                 I = c(V_A / mean^2)) %>%
      filter(is.na(var_sire) == F)
    
    REML.results2 <- REML.results.all %>%
      dplyr::select(c(V_A:I)) %>%
      pivot_longer(cols = c(1:4), names_to = "var_comp", values_to = "est")
    REML.results2$var_comp <- factor(REML.results2$var_comp, levels = c('V_P', 'V_A', 'h2', 'I')) #change order of VA and VP
    levels(REML.results2$var_comp) <- c('V[Phenotypic]', 'V[Additive]', 'h^2', 'I[A]') #change labels for VA and VP
    
    output = list(AllVarComps = REML.results.all, VarCompsLong = REML.results2)
    
    return(output)
  }
}


# Randomly sample n groups 
sample_n_groups <- function(grouped_df, size, replace = TRUE, weight=NULL) {  
  grp_var <- grouped_df %>% 
    groups %>%
    unlist %>% 
    as.character
  random_grp <- grouped_df %>% 
    summarise() %>% 
    sample_n(size, replace, weight) %>% 
    mutate(unique_id = 1:NROW(.))
  grouped_df %>% 
    right_join(random_grp, by = grp_var, relationship = 'many-to-many')
}


# Bootstrap all variance components
bootstrap_varcomps <- function(data, design, trait, n_iter = 1000){
  
  # trait mean
  trait_mean <- data %>% rename(Y = all_of(trait)) %>% summarize(mean = mean(Y, na.rm = T)) %>% as.numeric()
  
  n_groups = length(unique(data$sire))
  
  #list to store results
  datalist = vector("list", length = n_iter) 
  
  #Loop to run multiple times
  for (i in 1:n_iter) {
    # sample new data
    sampled <- data %>% group_by(sire) %>% sample_n_groups(size = n_groups) %>%
      rename(Y = {{ trait }})
    
    if(design == 'fs'){
      # fit mixed model on sampled data
      mixedsamp <- lmer( Y ~ (1 | unique_id/dam), data = sampled)
      
      VarComps2 <-   # calculate variance components
        var_from_model(model = mixedsamp, mean = trait_mean, design = design)
    }
    
    if(design == 'phs'){
      # fit mixed model on sampled data
      mixedsamp <- lmer( Y ~ (1 | unique_id), data = sampled)
      
      VarComps2 <-   # calculate variance components
        var_from_model(model = mixedsamp, mean = trait_mean, design = design)
    }
    
    datalist[[i]] <- VarComps2$AllVarComps
  }
  
  datalist
  
}


# Summarize outputs from bootstrapping
boot_summary <- function(bootout){
  #summarize results for VA
  summ_A <- bootout %>% summarize(  
    mean = mean(V_A),
    sd = sd(V_A),
    se = sd(V_A)/sqrt(n()),
    CI.low = quantile(V_A, 0.025),
    CI.high = quantile(V_A, 0.975)) %>%
    mutate(var_comp = c('VA'))
  
  #summarize results for VP
  summ_P <- bootout %>% summarize(  
    mean = mean(V_P),
    sd = sd(V_P),
    se = sd(V_P)/sqrt(n()),
    CI.low = quantile(V_P, 0.025),
    CI.high = quantile(V_P, 0.975)) %>%
    mutate(var_comp = c('VP'))
  
  #summarize results for h2
  summ_h <- bootout %>% summarize(  
    mean = mean(h2),
    sd = sd(h2),
    se = sd(h2)/sqrt(n()),
    CI.low = quantile(h2, 0.025),
    CI.high = quantile(h2, 0.975)) %>%
    mutate(var_comp = c('h2'))
  
  #summarize results for h2
  summ_I <- bootout %>% summarize(  
    mean = mean(I),
    sd = sd(I),
    se = sd(I)/sqrt(n()),
    CI.low = quantile(I, 0.025),
    CI.high = quantile(I, 0.975)) %>%
    mutate(var_comp = c('I'))
  
  rbind(summ_A,summ_P,summ_h, summ_I)  
}



## Days until Diapause ----


### Import data of diapause experiment ----
y.n = read_csv("data/diapause_herit_data.csv")

### Data Clean up

y.n <- y.n %>%
  filter(is.na(diapause) == F ) %>% # filter out non-trials
  select(1:7) # remove columns for each day
y.n$sire <- as.factor(y.n$sire)
str(y.n)

y.n %>% group_by(treatment, sire) %>% summarise(
  noff = n()) %>% summarise(
    range = range(noff),
    nsire = n(),
    mean = mean(noff)
  )

y.n %>% group_by(treatment) %>% summarise(
  noff = n()) # sample size

y.n %>% group_by(sire, dam) %>% summarise(
  noff = n()) #sample size of full-sib families

# sample size per sire
y.n %>% group_by(sire, treatment) %>% 
  filter(is.na(dead) == TRUE) %>%
  summarise(
    noff = n())

# sires with too few offspring to be included
exclude_small_families <- y.n %>% group_by(sire, treatment) %>% 
  filter(is.na(dead) == TRUE) %>%
  summarise(
    noff = n()) %>%
  filter(noff < 2)

short <- y.n %>% # create dataset for short treatment
  filter(treatment == "S" & is.na(dead) == TRUE) %>% # filter out dead beetles
  anti_join(., exclude_small_families, by = c("sire", "treatment")) # remove sires with too few offspring
hist(short$time, breaks = 50)

long <- y.n %>% # create dataset for long treatment
  filter(treatment == "L" & is.na(dead) == TRUE) %>% # filter out dead beetles
  anti_join(., exclude_small_families, by = c("sire", "treatment")) # remove sires with too few offspring
hist(long$time, breaks = 50)
str(long)

short %>% group_by(sire) %>% 
  summarise(
    noff = n()) %>%
  summarise(
    max = max(noff),
    min = min(noff),
    nsire = n(),
    mean = mean(noff)
  )
long %>% group_by(sire) %>% 
  summarise(
    noff = n()) %>%
  summarise(
    max = max(noff),
    min = min(noff),
    nsire = n(),
    mean = mean(noff)
  )
# total individuals measured per treatment
short %>% summarise(n())
long %>% summarise(n())

# # range of days to diapause in short and long
# y.n %>% group_by(treatment) %>%
#   summarise(min = min(time, na.rm = T),
#             max = max(time, na.rm = T),
#             mean = mean(time, na.rm = T))

### Random effects Models ----
#Long day treatment
days_meanL <- long %>% summarise(mean = mean(time)) %>% as.numeric()
mixedLong1 <- lmer(time ~ (1 | sire), data = long)
summary(mixedLong1)
# confint(mixedLong1, method = "boot", oldNames = F)
confint(mixedLong1, method = "profile", oldNames = F)
rand(mixedLong1)
days_res_L <- var_from_model(model = mixedLong1, mean = days_meanL, design = 'phs')
days_res_forplot_L <- days_res_L$VarCompsLong %>%
  mutate(treatment = "long")

#Short day treatment
days_meanS <- short %>% summarise(mean = mean(time)) %>% as.numeric()
mixedShort1 <- lmer(time ~ (1 | sire), data = short)
summary(mixedShort1)
# confint(mixedShort1, method = "boot")
confint(mixedShort1, method = "profile")
rand(mixedShort1)
days_res_S <- var_from_model(mixedShort1, days_meanS, design = 'phs')
days_res_forplot_S <- days_res_S$VarCompsLong %>%
  mutate(treatment = "short")

###Days Bootstrap confidence intervals ----
#long
# data, design, trait,
boot_days_L <- bootstrap_varcomps(data = long, design = 'phs', trait = 'time') %>% 
  bind_rows() %>% filter(is.na(var_sire) == F) #combine list into one dataframe, remove NAs
boot_sum_days_L <- boot_summary(boot_days_L) %>%
  mutate(treatment = "long")

#short
boot_days_S <- bootstrap_varcomps(data = short, design = 'phs', trait = 'time') %>% 
  bind_rows() %>% filter(is.na(var_sire) == F) #combine list into one dataframe, remove NAs
boot_sum_days_S <- boot_summary(boot_days_S) %>%
  mutate(treatment = "short")

###Plots of Results ----
days_res_comb <- rbind(days_res_forplot_L, days_res_forplot_S) # means from REML

days_boot_sum_both <- rbind(boot_sum_days_L,boot_sum_days_S) # confidence intervals from bootstrap
days_boot_sum_both$var_comp <- factor(days_boot_sum_both$var_comp, levels = c('VP', 'VA', 'h2', 'I')) #change order of VA and VP
levels(days_boot_sum_both$var_comp) <- c('V[Phenotypic]', 'V[Additive]', 'h^2', 'I[A]') #change labels for VA and VP

Vars <- ggplot(data = days_boot_sum_both %>% filter(var_comp != 'h^2' & var_comp != 'I[A]'), 
               aes(x = treatment, y = mean, color = treatment)) +
  #sd error bars
  # geom_linerange(aes(ymax = mean + sd, ymin = mean - sd),
  #                size = 1, lty = 'solid') +
  #conf int error bars
  geom_linerange(aes(ymax = CI.high, ymin = CI.low),
                 size = 1, lty = 'solid') +
  geom_point(data = days_res_comb %>% filter(var_comp != 'h^2' & var_comp != 'I[A]'), 
             aes(x = treatment, y = est, shape = treatment), fill = 'white',
             size = 3) +
  facet_wrap( ~ var_comp, labeller = label_parsed, scales = "fixed") +
  labs(y = "Variance", x = "Environment") +
  scale_x_discrete(labels = c('long' = "North\n(Home)", 'short' = "South\n(Away)"))+
  scale_shape_manual(labels = c('long' = "Home", 'short' = "Away"), 
                     values = c('long' = 16, 'short' = 21)) +
  scale_color_manual(labels = c('long' = "Home", 'short' = "Away"), 
                     values = c('long' = '#006b4d', 'short' = '#00eba9')) +
  theme_bw() +
  theme(text = element_text(size = 18), #axis.ticks.x = element_blank(),
        strip.placement = "outside", strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'none')
Vars

Herit <- ggplot(data = days_boot_sum_both %>% filter(var_comp == 'h^2'), 
                aes(x = treatment, y = mean, color = treatment)) +
  #sd error bars
  # geom_linerange(aes(ymax = mean + sd, ymin = mean - sd),
  #                size = 1, lty = 'solid') +
  # scale_y_continuous(breaks = c(0,0.5,1), limits = c(-0.05, 1.14)) +
  #conf int error bars
  geom_linerange(aes(ymax = CI.high, ymin = CI.low),
                 size = 1, lty = 'solid') +
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(-0.05, 1.49)) +
  #h2 means
  geom_point(data = days_res_comb %>% filter(var_comp == 'h^2'), 
             aes(x = treatment, y = est, shape = treatment), fill = 'white',
             size = 3) +
  #formatting
  facet_wrap( ~ var_comp, labeller = label_parsed, scales = "fixed") +
  labs(y = "Heritability", x = "Environment") +
  scale_x_discrete(labels = c('long' = "North\n(Home)", 'short' = "South\n(Away)")) +
  scale_shape_manual(values = c('long' = 16, 'short' = 21)) +
  scale_color_manual(labels = c('long' = "Home", 'short' = "Away"), 
                     values = c('long' = '#006b4d', 'short' = '#00eba9')) +
  theme_bw() +
  theme(text = element_text(size = 18), #axis.ticks.x = element_blank(),
        strip.placement = "outside", strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'none')
Herit

Evolv <- ggplot(data = days_boot_sum_both %>% filter(var_comp == 'I[A]'), 
                aes(x = treatment, y = mean, color = treatment)) +
  #sd error bars
  # geom_linerange(aes(ymax = mean + sd, ymin = mean - sd),
  #                size = 1, lty = 'solid') +
  #conf int error bars
  geom_linerange(aes(ymax = CI.high, ymin = CI.low),
                 size = 1, lty = 'solid') +
  geom_point(data = days_res_comb %>% filter(var_comp == 'I[A]'), 
             aes(x = treatment, y = est, shape = treatment), fill = 'white',
             size = 3) +
  facet_wrap( ~ var_comp, labeller = label_parsed, scales = "fixed") +
  labs(y = "Evolvability", x = "Environment") +
  scale_x_discrete(labels = c('long' = "North\n(Home)", 'short' = "South\n(Away)")) +
  scale_shape_manual(values = c('long' = 16, 'short' = 21)) +
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(-0.05, 1.14)) +
  scale_color_manual(labels = c('long' = "Home", 'short' = "Away"), 
                     values = c('long' = '#006b4d', 'short' = '#00eba9')) +
  theme_bw() +
  theme(text = element_text(size = 18), #axis.ticks.x = element_blank(),
        strip.placement = "outside", strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'none')
Evolv

ggarrange(Vars, Herit, Evolv, ncol = 3, widths = c(.8, 0.5, 0.5), labels = 'AUTO')

ggsave("plots/Fig5_diapause_h2.png", width = 9, height = 4, units = "in")



### Plot of variation among families ----

resim_short <- REsim(mixedShort1) %>% mutate(treat = 'short') %>% #short dataset
  mutate(addIntercept = mean + fixef(mixedShort1)) %>%
  mutate(conf = sd * qnorm(1-((1-0.95)/2)))

resim_long <- REsim(mixedLong1) %>% mutate(treat = 'long') %>% #long dataset
  mutate(addIntercept = mean + fixef(mixedLong1)) %>% 
  # arrange(addIntercept) %>%
  # mutate(groupID_order = row_number()) %>%
  # mutate(groupID2 = factor(groupID, levels = groupID)) %>%
  mutate(conf = sd * qnorm(1-((1-0.95)/2)))
plotREsim(resim_long, stat = 'mean')
str(resim_long)

order_key <- resim_long %>%
  arrange(addIntercept) %>%
  mutate(groupID_order = seq(1:nrow(resim_long))) %>%
  select(groupID, groupID_order)

resim_short <- resim_short %>%
  left_join(., order_key, by = 'groupID')
resim_long <- resim_long %>%
  left_join(., order_key, by = 'groupID') 
range(resim_long$addIntercept)
range(resim_short$addIntercept)


long.plot <- long %>% 
  left_join(., order_key, by = c('sire' = 'groupID'))
short.plot  <- short %>% 
  left_join(., order_key, by = c('sire' = 'groupID'))


varAmongFams <- ggplot(data = resim_long, aes(x = groupID_order, y = addIntercept)) +
  # geom_point(aes(fill = treat), size = 3, shape = 21) +  #For means, without sd
  #add background points
  geom_point(data = long.plot, 
             aes(x = groupID_order, y = time), 
             position = position_jitter(width = 0.15), 
             size = 1, shape = 19, alpha = 0.5, color = '#006b4d') +
  geom_point(data = short.plot, aes(x = groupID_order, y = time), 
             position = position_jitter(width = 0.15), 
             size = 1, shape = 19, alpha = 0.5, color = '#00eba9', fill = 'white') +
  
  #add family means and error bars
  # geom_pointrange(aes(ymin = addIntercept - sd, ymax = addIntercept + sd, shape = treat), #For means with sd
  #                 fatten = 2, size = 1, fill = 'white') + #shape = 21) +
  # geom_pointrange(data = resim_short, aes(x = groupID, ymin = addIntercept - sd, ymax = addIntercept + sd, shape = treat),
  # fatten = 2, size = 1, fill = 'white') + #, shape = 21) +
  
  #means with CI
  geom_pointrange(aes(x = groupID_order, ymin = addIntercept - conf, ymax = addIntercept + conf, 
                      shape = treat), #For means with CI
                  fatten = 4, size = 0.5, color = '#006b4d', fill = 'white', 
                  position = position_nudge(x = -0.1)) + #shape = 21) +
  geom_pointrange(data = resim_short, 
                  aes(x = groupID_order, ymin = addIntercept - conf, ymax = addIntercept + conf, 
                      shape = treat),
                  fatten = 4, size = 0.5, color = '#00eba9', fill = 'white', 
                  position = position_nudge(x = 0.1)) + #, shape = 21) +
  #formatting
  # lims(y = c(0,31)) +
  labs(x = "Rank order of families", y = "Days until diapause", shape = "Environment") +
  scale_shape_manual(labels = c('long' = "North (home)", 'short' = "South (away)"), 
                     values = c(16, 21))+# c('long' = , 'short' = '2')) +
  theme_classic() +
  theme(text = element_text(size = 18), axis.ticks.x = element_blank(),
        axis.text.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'bottom')
varAmongFams

ggsave("plots/Fig4_variation_families_new.png", width = 7.21, height = 4.17, units = "in")


