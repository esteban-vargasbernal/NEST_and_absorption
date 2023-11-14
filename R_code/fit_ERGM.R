library(ergm.ego)
library(parallel)
library(dplyr)
library(tidyr)  
library(knitr)
library(ggplot2)
library(reshape2)
library(gplots)

load("NESTegodata_both.RData") # this is stored as the egodata 'NEST' and contains the NEST data

# this is an egor file with distinct sample weights in egoWt
NEST_1 <- as_egor.egodata(NEST) 

# this is an egor file with same sample weights in eogWt
egos_temp <- NEST$egos
egos_temp$egoWt <- 1
NEST2 <- list(egos = egos_temp, alters = NEST$alters, egoWt = egos_temp$egoWt, egoIDcol = NEST$egoIDcol)
NEST_2 <- as_egor.egodata(NEST2) 

# We fit two ERGMs; one for the data with distinct sample weights and one for the data with no distinct sample weights 

N <- 3000 

model1 <- ergm.ego(NEST_1~                           
                     nodefactor("ethcat", levels=c("white","black")) +
                     nodefactor("agecat", levels = c("age.gp1","age.gp2","age.gp3")) +
                     nodefactor("treat", levels = c("Tr")) +
                     nodefactor("Trade_sex", levels = c("Trade")) + 
                     nodefactor("Group_sex", levels = c("Group")) + 
            
                     
                     nodematch("ethcat", diff = TRUE) +
                     nodematch("agecat", diff=TRUE) +
                     nodematch("treat", diff = TRUE)  + 
                     
                    degree(1, by = "ethcat") 
                   ,
                   
                   control=control.ergm.ego(ergm.control=control.ergm(
                     MCMLE.maxit=100,
                     MCMLE.trustregion=10000, 
                     MPLE.max.dyad.types=1e5,
                     MCMC.interval=2e5
                   ), ppopsize=N),
                   popsize=1)


model2 <- ergm.ego(NEST_2~                           
                     nodefactor("ethcat", levels=c("white","black")) +
                     nodefactor("agecat", levels = c("age.gp1","age.gp2","age.gp3")) +
                     nodefactor("treat", levels = c("Tr")) +
                     nodefactor("Trade_sex", levels = c("Trade")) + 
                     nodefactor("Group_sex", levels = c("Group")) + 
                     
                     
                     nodematch("ethcat", diff = TRUE) +
                     nodematch("agecat", diff=TRUE) +
                     nodematch("treat", diff = TRUE)  + 
                     
                     degree(1, by = "ethcat") 
                   ,
                   
                   control=control.ergm.ego(ergm.control=control.ergm(
                     MCMLE.maxit=100,
                     MCMLE.trustregion=10000, 
                     MPLE.max.dyad.types=1e5,
                     MCMC.interval=2e5
                   ), ppopsize=N),
                   popsize=1)

save(model1, file = "model_dist_weights.RData") 
save(model2, file = "model_same_weights.RData") 


# get degree summaries 

degreedist(model2$network,plot = TRUE, by = 'treat')
GOFD2 = gof(model2, GOF = 'degree')
barplot(colMeans(GOFD2$psim.deg))
barplot(degreedist(NEST_1))




## With this code, we take the attributes of all the nodes in the sampled networks from our ERGM


# change if needed 

load("model_dist_weights.RData")
model <- model1 
nsims <- 1



library(statnet.common)
formula <- nonsimp_update.formula(model$formula,model$newnetwork~edges+.,from.new=TRUE)
sim <- simulate(formula,coef = unname(coef(model)),nsim = nsims,control=control.simulate.formula(MCMC.interval = 1e6))

edges.ls <- vector(mode="list",length = nsims)

for(i in 1:nsims){
  len.tmp <- length(sim[[i]]$mel)
  edges.tmp <- matrix(NA, nrow =  2*len.tmp, ncol = 2)
  for(j in 1:len.tmp){
    edges.tmp[j,] <- c(sim[[i]]$mel[[j]]$outl, sim[[i]]$mel[[j]]$inl)
    edges.tmp[len.tmp+j,] <- c(sim[[i]]$mel[[j]]$inl, sim[[i]]$mel[[j]]$outl)
  }
  edges.ls[[i]] <- edges.tmp
}


for (i in 1:nsims) {
  edge.tmp <- edges.ls[[i]]
  write.csv(edge.tmp,paste0("sims/edges/edges_sim_",as.character(i),".csv"))
}



N <- sim[[1]]$gal$n
deltav <- matrix(0,nrow = N, ncol = 1)

for(i in 1:N){
  deltav[i,1] <- ifelse(sim[[1]]$val[[i]]$treat=="Tr",1,1/5)
}

write.csv(deltav,paste0("sims/attributes/delta_vector.csv"))



sim_i = 1

N <- length(sim[[sim_i]][["val"]])

summary_df_0 <- data.frame("node" = 1:N, "ethcat" = rep(0,N), "agecat" = rep(0,N), "treat" = rep(0,N), "Trade_sex" = rep(0,N),
                           "Group_sex" = rep(0,N), "degree" = rep(0,N),
                           "age.gp1" = rep(0,N), "age.gp2" = rep(0,N), "age.gp3" = rep(0,N), "age.gp4" = rep(0,N),
                           "black" = rep(0,N), "white" = rep(0,N), "other" = rep(0,N),
                           "Tr" = rep(0,N), "NoTr" = rep(0,N), "Trade" = rep(0,N), "NoTrade" = rep(0,N),
                           "Group" = rep(0,N), "NoGroup" = rep(0,N))


for(i in 1:N){
  summary_df_0[i,"ethcat"] <- sim[[sim_i]][["val"]][[i]]$ethcat
  summary_df_0[i, "agecat"] <- sim[[sim_i]][["val"]][[i]]$agecat
  summary_df_0[i, "treat"] <- sim[[sim_i]][["val"]][[i]]$treat
  summary_df_0[i, "Trade_sex"] <- sim[[sim_i]][["val"]][[i]]$Trade_sex
  summary_df_0[i, "Group_sex"] <- sim[[sim_i]][["val"]][[i]]$Group_sex
}


for (sim_i in 1:nsims) {
  
  Ne <- length(sim[[sim_i]][["mel"]])
  edges.tmp <- matrix(NA, nrow = 2*Ne, ncol = 2)
  for(j in 1:Ne){
    edges.tmp[j,] <- c(sim[[sim_i]][["mel"]][[j]]$outl, sim[[sim_i]][["mel"]][[j]]$inl)
    edges.tmp[Ne+j,] <- c(sim[[sim_i]][["mel"]][[j]]$inl, sim[[sim_i]][["mel"]][[j]]$outl)
  }
  
  summary_df <- summary_df_0
  for(i in 1:N){
    for(e in 1:(2*Ne)){
      if(edges.tmp[e,1]==i){
        j = edges.tmp[e,2]
        summary_df[i,"degree"] = summary_df[i,"degree"] + 1
        if(sim[[sim_i]][["val"]][[j]][["ethcat"]] == "black"){summary_df[i,"black"] = summary_df[i,"black"]+1}
        if(sim[[sim_i]][["val"]][[j]][["ethcat"]] == "white"){summary_df[i,"white"] = summary_df[i,"white"]+1}
        if(sim[[sim_i]][["val"]][[j]][["ethcat"]] == "other"){summary_df[i,"other"] = summary_df[i,"other"]+1}
        if(sim[[sim_i]][["val"]][[j]][["agecat"]] == "age.gp1"){summary_df[i,"age.gp1"] = summary_df[i,"age.gp1"]+1}
        if(sim[[sim_i]][["val"]][[j]][["agecat"]] == "age.gp2"){summary_df[i,"age.gp2"] = summary_df[i,"age.gp2"]+1}
        if(sim[[sim_i]][["val"]][[j]][["agecat"]] == "age.gp3"){summary_df[i,"age.gp3"] = summary_df[i,"age.gp3"]+1}
        if(sim[[sim_i]][["val"]][[j]][["agecat"]] == "age.gp4"){summary_df[i,"age.gp4"] = summary_df[i,"age.gp4"]+1}
        if(sim[[sim_i]][["val"]][[j]][["treat"]] == "Tr"){summary_df[i,"Tr"] = summary_df[i,"Tr"]+1}
        if(sim[[sim_i]][["val"]][[j]][["treat"]] == "NoTr"){summary_df[i,"NoTr"] = summary_df[i,"NoTr"]+1}
        if(sim[[sim_i]][["val"]][[j]][["Trade_sex"]] == "Trade"){summary_df[i,"Trade"] = summary_df[i,"Trade"]+1}
        if(sim[[sim_i]][["val"]][[j]][["Trade_sex"]] == "NoTrade"){summary_df[i,"NoTrade"] = summary_df[i,"NoTrade"]+1}
        if(sim[[sim_i]][["val"]][[j]][["Group_sex"]] == "Group"){summary_df[i,"Group"] = summary_df[i,"Group"]+1}
        if(sim[[sim_i]][["val"]][[j]][["Group_sex"]] == "NoGroup"){summary_df[i,"NoGroup"] = summary_df[i,"NoGroup"]+1}
      }
    }
  }
  
  write.csv(summary_df, paste0("sims/attributes/summary_sim",as.character(sim_i),".csv"))
  
  desc_nodes <- summary_df %>% group_by(degree) %>% count() 
  desc_nodes <- desc_nodes[2:length(desc_nodes$degree), ] 
  Nt <- sum(desc_nodes$n)
  desc_nodes <- desc_nodes %>% mutate(Percent = round(n/Nt,5))
  
  write.csv(desc_nodes,paste0("sims/attributes/degree_dist_sim_",as.character(sim_i),".csv"))
  
}



