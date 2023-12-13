
DD_data <- read_excel('DD_model.xlsx', sheet = 'combined')


ggplot(DD_data, aes (x = dayofyear, y = CUMDD11, color = region)) +
  geom_line() +
  geom_vline(xintercept = 209, color = 'red') +
  geom_vline(xintercept = 250, color = 'blue') +
  theme_bw()

##Long treatment - Delta 
daystodiapause_long <- meansDiapausingPops_df %>% mutate(core_edge = c(rep(c(rep("core", 4), rep('edge', 4)),2))) %>% filter(treatment == 'long') %>% filter(is.na(rate) == F) %>%
  mutate(dayofdiapause = rate + 209) %>% arrange(dayofdiapause)


ddnorth <- ggplot(DD_data %>% filter(region == 'delta'), aes (x = dayofyear, y = CUMDD11)) +
  geom_line() +
  geom_vline(xintercept = 209, color = 'red', size = 1.5) + #day of year when treatment occurs
  geom_vline(xintercept = 287, color = 'blue', size = 1.5) +  # day of first frost
  geom_vline(data = daystodiapause_long, aes(xintercept = dayofdiapause, color = core_edge), size = 1) +
  labs(x = "Day of Year", y = "Cumulative Degree Days", color = "Region", title = "Degree days in Delta (north)") +
  theme_bw()
ddnorth

##Short treatment - Topock Marsh
daystodiapause_short <- meansDiapausingPops_df %>% mutate(core_edge = c(rep(c(rep("core", 4), rep('edge', 4)),2))) %>% filter(treatment == 'short') %>% filter(is.na(rate) == F) %>%
  mutate(dayofdiapause = rate + 250) %>% arrange(dayofdiapause)


ddsouth <- ggplot(DD_data %>% filter(region == 'topock'), aes (x = dayofyear, y = CUMDD11)) +
  geom_line() +
  geom_vline(xintercept = 250, color = 'red', size = 1.5) + #day of year when treatment occurs
  geom_vline(xintercept = 365, color = 'blue', size = 1.5) +  # day of first frost
  geom_vline(data = daystodiapause_short, aes(xintercept = dayofdiapause, color = core_edge), size = 1) +
  labs(x = "Day of Year", y = "Cumulative Degree Days", color = "Region", title = "Degree days in Topock Marsh (south)") +
  theme_bw()
ddsouth

ggarrange(ddnorth, ddsouth, ncol = 2, common.legend = T, legend = 'bottom')
