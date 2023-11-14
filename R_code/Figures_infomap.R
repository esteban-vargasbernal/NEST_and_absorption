library(tidyverse)
library(ggplot2)
library(dplyr)


###### sizes of largest PrEP and no PrEP communities

node_comm = read.csv('Summaries/summary_all_node_comm.csv')
node_comm_tmp <- node_comm %>% filter(type_comm_no_tr=='large') %>% group_by(type_comm_tr,net,treat) %>% count() 

n_lccs <- node_comm %>% filter(type_comm_tr == 'large',type_comm_no_tr == 'large') %>% 
  group_by(net) %>% summarise(large_PrEP = max(size_comm_tr), large_no_PrEP = max(size_comm_no_tr))
n_lccs <- data_frame(sizes = c(n_lccs$large_no_PrEP, n_lccs$large_PrEP), Scenario = c(rep('no-PrEP',nsim),c(rep('PrEP',nsim))))

g_lccs <- n_lccs %>%
  ggplot(aes(x = sizes, fill = Scenario)) +
  geom_histogram(position = 'dodge', colour = 'black', bins = 20) +
  labs(x = 'Sizes of large communities') +
  labs(y = 'Percent of networks') +
  theme_bw() + 
  scale_fill_manual(values = c('grey40','grey')) + 
  theme(text = element_text(size = 15))
theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))

g_lccs

ggsave(g_lccs, file = paste0('Figures/largest_community_sizes.eps'),device = 'eps')

n_lccs %>% group_by(Scenario) %>% summarise(mean = mean(sizes),sd = sd(sizes))


###### Number of small PrEP communities inside large no-PrEP community 

str <- read.csv('Summaries/summary_all_comm_in_large_no_PrEP_comm.csv')

n_small_comm <- str %>% filter(comm > 1) %>% group_by(net) %>% count(net)

g_small_comm <- n_small_comm %>% 
  ggplot(aes(x = n)) + 
  geom_histogram(fill = 'grey', colour = 'black') + 
  labs(x = 'Number of small PrEP communities \n inside the large no-PrEP community') + 
  labs(y = 'Percent of networks') + 
  theme_bw() +
  theme(text = element_text(size = 15))

g_small_comm

ggsave(g_small_comm, file = paste0('Figures/number_small_PrEP_comm.eps'),device = 'eps')

n_small_comm %>% ungroup() %>% summarise_at(vars(n), list(mean=mean, sd =sd))



##### Number of PrEP nodes in large no-PrEP community that belong to small PrEP communities 

prop_PrEP <- node_comm_tmp %>% filter(treat == 'Tr') %>% group_by(net) %>% mutate(prop = n/sum(n))

g_prop_PrEP <- prop_PrEP %>% filter(type_comm_tr == 'small') %>%
  ggplot(aes(x = n)) + 
  geom_histogram(position = 'dodge', fill = "grey", colour = 'black') +
  labs(x = 'Number of PrEP nodes of the \n large no-PrEP comm that belong to small PrEP comm') + 
  labs(y = 'Percent of networks') +
  theme_bw() + 
  scale_fill_manual(values = c("grey","grey")) +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))

g_prop_PrEP

ggsave(g_prop_PrEP, file = paste0('Figures/n_PrEP_nodes_in_small_PrEP_comm_vs_large_comm.eps'),device = 'eps')

prop_PrEP %>% filter(type_comm_tr == 'small') %>% ungroup() %>% 
  summarise_at(vars(n),list(mean = mean, sd = sd))

##### Number of PrEP nodes in large no-PrEP community that belong to the large PrEP community 

prop_PrEP <- node_comm_tmp %>% filter(treat == 'Tr') %>% group_by(net) %>% mutate(prop = n/sum(n))

g_prop_PrEP <- prop_PrEP %>% filter(type_comm_tr == 'large') %>%
  ggplot(aes(x = n)) + 
  geom_histogram(position = 'dodge', fill = "grey", colour = 'black') +
  labs(x = 'Number of PrEP nodes of the \n large no-PrEP comm that belong to small PrEP comm') + 
  labs(y = 'Percent of networks') +
  theme_bw() + 
  scale_fill_manual(values = c("grey","grey")) +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))

g_prop_PrEP

ggsave(g_prop_PrEP, file = paste0('Figures/n_PrEP_nodes_in_large_PrEP_comm_vs_large_comm.eps'),device = 'eps')

prop_PrEP %>% filter(type_comm_tr == 'large') %>% ungroup() %>% 
  summarise_at(vars(n),list(mean = mean, sd = sd))


####################### SMALL COMMUNITY SIZES 

cluster <- read.csv('Summaries/summary_all_cluster.csv')

g_small_comm <- cluster %>% group_by(net,comm,scenario) %>% count() %>% filter(n<1000) %>%  #group_by(n, scenario,net) %>% count(n) %>%  mutate(prop = nn/sum(nn)) #%>% #%>% group_by(n) %>% count() %>%
  ggplot(aes(x = n, fill = scenario)) +
  geom_bar(aes(y = 100*..prop..),position = 'dodge') +
  labs(x = 'community size', y = 'Percent of small communities') +
  ggtitle(paste0('Sizes of small communities')) + 
  scale_x_discrete(limits = c(1:9)) +
  theme_bw() +
  scale_fill_manual(name = "Scenario", labels =c("no_PrEP", "PrEP"),values = c("grey40","grey")) +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))

g_small_comm

ggsave(g_small_comm, file = 'Figures/sizes_of_small_comm.eps',device = 'eps')


prop_cluster_small <- cluster %>% group_by(net,comm,scenario) %>% count() %>% filter(n<1000)  %>%  group_by(n,scenario) %>% count() %>% group_by(scenario) %>% mutate(prop = 100* nn/sum(nn))   
prop_cluster_small_by_net <- cluster %>% group_by(net,comm,scenario) %>% count() %>% filter(n<1000)  %>%  group_by(n,scenario,net) %>% count() %>% group_by(net,scenario) %>% mutate(prop = nn/sum(nn))   
cluster %>% group_by(net,comm,scenario) %>% count() %>% filter(n<1000)  %>%  group_by(n,scenario,net) %>% count() %>% group_by(n,scenario) %>% summarise(mean = mean(nn))   


############# NEIGHBORHOOD OF RADIUS 4

ego4 = read.csv('Summaries/summary_all_ego4.csv')

g_ego4_prop <- ego4 %>% 
  ggplot(aes(x=prop_trt, fill = type_PrEP_comm)) +
  geom_histogram(aes(y =c(..density..[..group..==1],
                          ..density..[..group..==2])),
                 position = 'dodge',bins = 20, colour = 'black') +
  labs(x = 'Proportion of PrEP nodes') +
  labs(y = 'Percent') + 
  theme_bw() +
  scale_fill_manual(name = "PrEP community", labels =c("Large", "Small"),values = c("grey40","grey")) +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))

g_ego4_prop

ggsave(g_ego4_prop, file = 'Figures/PrEP_neighbors.eps',device = 'eps')


ego4 %>%  group_by(type_PrEP_comm) %>% summarise(mean = mean(prop_trt), sd = sd(prop_trt))

######### display some ERGMs

i0 <-1
g_ego4_prop <- ego4 %>% filter(net==i0) %>%
  ggplot(aes(x=prop_trt, fill = type_PrEP_comm)) +
  geom_histogram(aes(y =c(..density..[..group..==1],
                          ..density..[..group..==2])),
                 position = 'dodge',bins = 20, colour = 'black') +
  labs(x = 'Proportion of PrEP nodes') +
  labs(y = 'Percent') + 
  theme_bw() +
  scale_fill_manual(name = "PrEP community", labels =c("Large", "Small"),values = c("grey40","grey")) +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))

g_ego4_prop

ggsave(g_ego4_prop, file = paste0('Figures/PrEP_neighbors_net',i0,'.eps'),device = 'eps')


ego4 %>% filter(net== i0) %>% group_by(type_PrEP_comm) %>% summarise(mean = mean(prop_trt), sd = sd(prop_trt))

######## degree inside the large no-PrEP comm


g_ego4_dg <- ego4  %>%
  ggplot(aes(x = degree, fill = type_PrEP_comm)) + 
  geom_histogram(aes(y = c(..density..[..group..==1],
                           ..density..[..group..==2])), position = 'dodge', bins = 20, colour = 'black') +
  labs(x = 'Degree') + 
  labs(y = 'Percent') + 
  theme_bw() +
  scale_fill_manual(name = "PrEP community", labels =c("Large", "Small"),values = c("grey40","grey")) +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))

g_ego4_dg

ggsave(g_ego4_dg, file = paste0('Figures/degree_by_PrEP_comm.eps'),device = 'eps')


ego4  %>% group_by(type_PrEP_comm) %>%
  summarise_at(vars(degree), list(mean=mean, sd = sd))

######## display some ERGMs

i0 = 1
g_ego4_dg <- ego4  %>% filter(net==1) %>%
  ggplot(aes(x = degree, fill = type_PrEP_comm)) + 
  geom_histogram(aes(y = c(..density..[..group..==1],
                           ..density..[..group..==2])), position = 'dodge', bins = 20, colour = 'black') +
  labs(x = 'Degree') + 
  labs(y = 'Percent') + 
  theme_bw() +
  scale_fill_manual(name = "PrEP community", labels =c("Large", "Small"),values = c("grey40","grey")) +
  theme(text = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))

g_ego4_dg

ggsave(g_ego4_dg, file = paste0('Figures/degree_by_PrEP_comm_net',i0,'.eps'),device = 'eps')


ego4  %>% filter(net==1) %>% group_by(type_PrEP_comm) %>%
  summarise_at(vars(degree), list(mean=mean, sd = sd))

