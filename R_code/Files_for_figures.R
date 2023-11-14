library(dplyr)
library(tidyr)
library(knitr)
library(ergm.ego)
library(ggplot2)

# chnage if needed
nsim <- 2

#########################
# Getting all partitions yielded by InfoMap

for(i0 in 1:nsim){
cluster_tr = read.csv(paste0('../Python_code/infomap/PrEP/partition_PrEP_sim_',i0,'.csv'), header = FALSE )
cluster_no_tr = read.csv(paste0('../Python_code/infomap/no_PrEP/partition_no_PrEP_sim_',i0,'.csv'), header= FALSE)

colnames(cluster_tr) <- c('comm')
colnames(cluster_no_tr) <- c('comm')

cluster_tr <- cluster_tr %>% mutate(scenario = 'tr') %>% mutate(node = 1:nrow(cluster_tr)) %>% mutate(net = rep(i0,nrow(cluster_tr)))
cluster_no_tr <- cluster_no_tr %>% mutate(scenario = 'no_tr') %>% mutate(node = 1:nrow(cluster_no_tr)) %>% mutate(net = rep(i0, nrow(cluster_tr)))

cluster_tmp <- union_all(cluster_tr, cluster_no_tr)

if(i0 == 1){cluster <- cluster_tmp} else{cluster <- union_all(cluster, cluster_tmp)}
}

write.csv(cluster,'Summaries/summary_all_cluster.csv')


################# summary all community

for(i0 in 1:nsim){
  node_comm_tmp <- read.csv(paste0('../Python_code/infomap/attributes/nodes_comm_sim_',i0,'.csv'))
  node_comm_tmp <- node_comm_tmp %>% mutate(net = i0)
  if(i0==1){node_comm <- node_comm_tmp}else{node_comm <- union_all(node_comm,node_comm_tmp)}
}

write.csv(node_comm,'Summaries/summary_all_node_comm.csv')

node_comm_tmp <- node_comm %>% filter(type_comm_no_tr=='large') %>% group_by(type_comm_tr,net,treat) %>% count() 

prop_PrEP <- node_comm_tmp %>% filter(treat == 'Tr') %>% group_by(net) %>% mutate(prop = n/sum(n))
prop_small <- node_comm_tmp %>% group_by(type_comm_tr,net) %>% summarise(n = sum(n)) %>% group_by(net) %>% mutate(prop = n/sum(n))
prop_PrEP_small <- node_comm_tmp %>% filter(type_comm_tr == 'small') %>% group_by(net) %>% mutate(prop = n/sum(n))

########################
# PrEP neighbors 

for(i0 in 1:nsim){
ego4_large = read.csv(paste0('../Python_code/infomap/attributes/PrEP_neighbors_large_comm_sim_',i0,'.csv'))
ego4_large <- ego4_large %>% mutate(type_PrEP_comm = 'large') %>% mutate(net = i0)
ego4_small = read.csv(paste0('../Python_code/infomap/attributes/PrEP_neighbors_small_comm_sim_',i0,'.csv'))
ego4_small <- ego4_small %>% mutate(type_PrEP_comm = 'small') %>% mutate(net = i0)
ego4_tmp <- union_all(ego4_small,ego4_large)
if(i0==1){ego4 <- ego4_tmp}else{ego4 <- union_all(ego4,ego4_tmp)}
}


write.csv(ego4, 'Summaries/summary_all_ego4.csv')

##########################
library(scales)

###############################################################


for(i0 in 1:nsim){
str_tmp = read.csv(paste0('../Python_code/infomap/attributes/PrEP_prop_per_comm_inside_large_no_PrEP_sim',i0,'.csv'))
str_tmp <- str_tmp %>% mutate(net = i0)
if(i0==1){str <- str_tmp}else{str = union_all(str, str_tmp)}
}

write.csv(str,'Summaries/summary_all_comm_in_large_no_PrEP_comm.csv')

# get proportion of PrEP nodes of large no-PrEP comm in the large PrEP comm (about 10%)

n_small_comm <- str %>% filter(comm > 1) %>% group_by(net) %>% count(net)


################################

# Summary of all nets

for(i0 in 1:nsim){
  summary_tmp <- read.csv(paste0('../Python_code/sims/attributes/summary_lcc_sim',i0,'.csv'))
  summary_tmp <- summary_tmp %>% mutate(net = rep(i0,nrow(summary_tmp)))
  if(i0==1){summary_lcc <-summary_tmp}else{summary_lcc <- union_all(summary_lcc, summary_tmp)}
}

write.csv(summary_lcc, 'Summaries/summary_all_lcc.csv')


################################################################

# SHORT RUNS STARTING AT LARGE NO-PREP COMMUNITY - READ-WRITE FILES

for(i0 in 1:nsim){
output_all <- read.csv('../Python_code/dynamics/summary_all_short_runs.csv')
node_comm = read.csv(paste0('../Python_code/infomap/attributes/nodes_comm_sim_',i0,'.csv'))

ICs_fixed <- read.csv(paste0('../Python_code/dynamics/initial_conditions/initial_conditions_sim_',i0,'.csv'))
node_tmp <-node_comm[node_comm$node %in% ICs_fixed$node,]
E0_comm2 <- output_all %>% filter(network==i0) %>% group_by(E0) %>% summarise(n_large = 2*sum(outbreak_type=='large')) 
colnames(E0_comm2) <- c('node','n_large') ## n_large is the number of long outbreaks among short runs
E0_comm <- inner_join(E0_comm2, node_tmp,by = 'node') %>% mutate(net = rep(i0, nrow(E0_comm2)))
if(i0==1){E0_comm_reinf <- E0_comm}else{E0_comm_reinf <- union_all(E0_comm_reinf,E0_comm)}
}

write.csv(E0_comm_reinf, 'Summaries/summary_all_large_outbreaks_short_runs.csv')

for(i0 in 1:nsim){ for(case in c('case_A','case_B')){
reinf_tmp <- read.csv(paste0('../Python_code/dynamics/long_runs/num_reinf_',case,'_sim_',i0,'.csv'))  
reinf_tmp <- reinf_tmp %>% mutate(net = rep(i0,nrow(reinf_tmp))) %>% mutate(case = case) 
if(i0==1 && case == 'case_A'){reinf <- reinf_tmp}else{reinf <- union_all(reinf, reinf_tmp)}
  }
}

write.csv(reinf, 'Summaries/summary_all_reinfections.csv')
