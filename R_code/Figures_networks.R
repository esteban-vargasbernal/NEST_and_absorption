library(tidyverse)
library(ggplot2)
library(dplyr)


#######################  DEGREE

summary_ergm_all <- read.csv('sims/attributes/summary_egocentric.csv')

##### Degree of egocentric sample

g_degree_ergm <- summary_ergm_all %>%  
  ggplot(aes(x= degree, fill = factor(treat))) + 
  geom_histogram(aes(y =100*c(..density..[..group..==1],
                              ..density..[..group..==2])),
                 position = 'dodge', bins = 10, colour = "black") + 
  labs(fill='Node type')+
  labs(x = 'Degree')+
  labs(y = 'Percent')+
  theme_bw() +
  scale_fill_manual(name = "Node type", labels = c("non-PrEP", "PrEP"),values = c("grey40","grey")) +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
g_degree_ergm

ggsave(g_degree_ergm, file = 'Figures/degree_egocentric.eps',device = 'eps')


summary_ergm_all  %>% group_by(treat) %>%
  summarise_at(vars(degree), list(mean = mean, sd = sd))

summary_ergm_all  %>%
  summarise_at(vars(degree), list(mean = mean, sd=sd))


#### largest connected component sizes 

node_comm = read.csv('Summaries/summary_all_node_comm.csv')
node_comm_tmp <- node_comm %>% filter(type_comm_no_tr=='large') %>% group_by(type_comm_tr,net,treat) %>% count() 
lcc_sizes <- node_comm %>% group_by(net) %>% count()

g_lcc_sizes <- lcc_sizes %>% 
  ggplot(aes(x = n)) + 
  geom_histogram(position = 'dodge', fill = "grey", colour = 'black') +
  labs(x = 'Largest component sizes') + 
  labs(y = 'Percent of networks') +
  theme_bw() + 
  scale_fill_manual(values = c("grey","grey")) +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))

g_lcc_sizes

ggsave(g_lcc_sizes, file = paste0('Figures/largest_component_sizes.eps'),device = 'eps')

lcc_sizes %>% ungroup() %>% summarise_at(vars(n),list(mean = mean, sd = sd))


######## degree of data by PrEP usage of all 100 networks

summary_data_all <- read.csv('Summaries/summary_all_lcc.csv')

g_degree_data <- summary_data_all %>% #filter(net == i0) %>% 
  ggplot(aes(x= degree, fill = treat)) + 
  geom_histogram(aes(y =100*c(..density..[..group..==1],
                              ..density..[..group..==2])),
                 position = 'dodge', bins = 10, colour = "black") + 
  labs(x = 'Degree')+
  labs(y = 'Percent')+
  theme_bw() +
  scale_fill_manual(name = "Node type", labels = c("non-PrEP", "PrEP"),values = c("grey40","grey")) +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))

g_degree_data

ggsave(g_degree_data, file = 'Figures/degree_data.eps',device = 'eps')


summary_data_all %>% group_by(treat) %>%
  summarise_at(vars(degree), list(mean = mean,sd=sd))

summary_data_all  %>% 
  summarise_at(vars(degree), list(mean = mean,sd=sd))


######## degree of data by PrEP usage of one network

i0 <- 1

g_degree_data <- summary_data_all %>% filter(net == i0) %>% 
  ggplot(aes(x= degree, fill = treat)) + 
  geom_histogram(aes(y =100*c(..density..[..group..==1],
                              ..density..[..group..==2])),
                 position = 'dodge', bins = 10, colour = "black") + 
  labs(x = 'Degree')+
  labs(y = 'Percent')+
  theme_bw() +
  scale_fill_manual(name = "Node type", labels = c("non-PrEP", "PrEP"),values = c("grey40","grey")) +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))

g_degree_data

ggsave(g_degree_data, file = 'Figures/degree_data_one_net.eps',device = 'eps')


summary_data_all %>% filter(net ==i0) %>% group_by(treat) %>%
  summarise_at(vars(degree), list(mean = mean,sd=sd))

summary_data_all  %>% filter(net==i0) %>%
  summarise_at(vars(degree), list(mean = mean,sd=sd))


