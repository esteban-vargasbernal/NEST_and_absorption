library(tidyverse)
library(ggplot2)
library(dplyr)


###### degrees of the initial conditions

output_all = read.csv(paste0('../Python_code/dynamics/summary_all_short_runs.csv'))


g_str_degree <- output_all %>%  group_by(network,E0,PrEP_comm_type) %>% summarise(degree = max(degree)) %>% ungroup(PrEP_comm_type) %>%
  ggplot(aes(y =degree, fill = factor(PrEP_comm_type))) +
  geom_boxplot() + 
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y = 'Degree') + 
  scale_fill_manual(name = 'Type of PrEP \n community \n of the PrEP IC', labels =c("Large", "Small"),values = c("grey60","white")) +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.y = element_text(size = 15))

g_str_degree

ggsave(g_str_degree, file = paste0('Figures/degree_initial_conditions.eps'),device = 'eps')

output_all %>%  group_by(network,E0,PrEP_comm_type) %>% summarise(degree = max(degree)) %>% ungroup()  %>% group_by(PrEP_comm_type) %>%
  summarise_at(vars(degree), list(mean = mean, sd = sd, median = median, max = max, min = min ))

###### Display 3 ERGMs

g_str_degree <- output_all %>% filter(network %in% 1:3) %>% group_by(network,E0,PrEP_comm_type) %>% summarise(degree = max(degree)) %>% ungroup(PrEP_comm_type) %>%
  ggplot(aes(x = factor(network), y =degree, fill = factor(PrEP_comm_type))) +
  geom_boxplot() + 
  labs(y = 'Degree') + 
  scale_fill_manual(name = 'Type of PrEP \n community \n of the PrEP IC', labels =c("Large", "Small"),values = c("grey60","white")) +
  theme_bw() +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))

g_str_degree

ggsave(g_str_degree, file = paste0('Figures/degree_initial_conditions_show_nets.eps'),device = 'eps')



##### Boxplot of long outbreaks by type of PrEP community

output_all = read.csv(paste0('../Python_code/dynamics/summary_all_short_runs.csv'))

g_str_disease <- output_all %>%  group_by(network,E0,outbreak_type,PrEP_comm_type) %>% count() %>% ungroup(PrEP_comm_type) %>% filter(outbreak_type=='large') %>%
  ggplot(aes(y =2*n, fill = factor(PrEP_comm_type) )) +
  geom_boxplot() + 
  labs(y = 'Percent of long outbreaks') + 
  #ggtitle('PrEP ICs inside the largest non-PrEP community') +
  scale_fill_manual(name = 'Type of PrEP \n community \n of the PrEP IC', labels =c("Large", "Small"),values = c("grey60","white")) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.y = element_text(size = 15))

g_str_disease

ggsave(g_str_disease, file = paste0('Figures/long_outbreaks_in_short_runs.eps'),device = 'eps')


output_all %>% group_by(network,E0,outbreak_type,PrEP_comm_type) %>% count() %>% ungroup()  %>% filter(outbreak_type=='large') %>% group_by(PrEP_comm_type) %>%
  summarise_at(vars(n), list(mean = mean, sd = sd, median = median, max = max, min = min ))

### Display 3 ERGMs

g_str_disease <- output_all %>% filter(network %in% 1:3) %>%  group_by(network,E0,outbreak_type,PrEP_comm_type) %>% count() %>% ungroup(PrEP_comm_type) %>% filter(outbreak_type=='large') %>%
  ggplot(aes(x = factor(network), y =2*n, fill = factor(PrEP_comm_type) )) +
  geom_boxplot() + 
  labs(y = 'Percent of long outbreaks') + 
  labs(x = 'Networks') + 
  #ggtitle('PrEP ICs inside the largest non-PrEP community') +
  scale_fill_manual(name = 'Type of PrEP \n community \n of the PrEP IC', labels =c("Large", "Small"),values = c("grey60","white")) +
  theme_bw() +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))

g_str_disease

ggsave(g_str_disease, file = paste0('Figures/long_outbreaks_in_short_runs_show_nets.eps'),device = 'eps')



####################### EQUILIBRIA FOR LONG OUTBREAKS

output_all = read.csv(paste0('../Python_code/dynamics/long_runs/summary_all_long_runs.csv')) ## generated in python


g_equi <- output_all  %>% filter(equilibria>0) %>% 
  ggplot(aes(x = case, y = equilibria, fill = factor(case))) + 
  geom_boxplot() + 
  labs(x = "case") +
  labs(y = "Prevalence pseudo-equilibrium") +
  theme_bw() +
  scale_fill_manual(name = "Case", labels =c('Case A', 'Case B'),values = c("grey40","grey","white")) +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))

g_equi

ggsave(g_equi, file = 'Figures/equilibria_long_runs.eps',device = 'eps')


### mean and sd
output_all %>%  group_by(case) %>%
  summarise_at(vars(equilibria), list(mean= mean, sd=sd))





########################### REINFECTIONS

num_reinf = read.csv(paste0('Summaries/summary_all_reinfections.csv'))

g_reinf <- num_reinf %>% filter(type_no_tr=='large', treat == 'PrEP') %>%
  ggplot(aes(x = case, y = n, fill = type_tr)) +
  geom_boxplot() + 
  labs(y = "Number of reinfections") + 
  theme_bw() +
  scale_fill_manual(name = " PrEP \n community", labels =c("Large", "Small"),values = c("grey40","grey")) +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.y = element_text(size = 15))


g_reinf

ggsave(g_reinf, file = paste0('Figures/reinfections.eps'),device = 'eps')

num_reinf %>% group_by(treat,case) %>% 
  summarise_at(vars(n), list(mean = mean, sd = sd, median = median, max = max, min = min ))

num_reinf %>% group_by(type_tr, treat) %>% 
  summarise_at(vars(n), list(mean = mean, sd = sd, median = median, min = min, max = max)) 



####### Display some ERGMs


g_reinf <- num_reinf %>% filter(type_no_tr=='large', treat == 'PrEP', net %in% 1:3, case == 'case_B') %>%
  ggplot(aes(x = factor(net), y = n, fill = type_tr)) +
  geom_boxplot() + 
  labs(y = "Number of reinfections") + 
  labs(x= 'Network') +
  ggtitle('Case B') +
  theme_bw() +
  scale_fill_manual(name = " PrEP \n community", labels =c("Large", "Small"),values = c("grey40","grey")) +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.y = element_text(size = 15),plot.title = element_text(hjust = 0.5))


g_reinf

ggsave(g_reinf, file = paste0('Figures/reinfections_show_nets.eps'),device = 'eps')




