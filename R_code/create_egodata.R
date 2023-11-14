library(dplyr)
library(tidyr)
library(knitr)
library(ergm.ego)
library(ggplot2)

# Reading the NEST data

table1.2019 <- read.csv("data/table1_2019.csv", header = TRUE)
table2.2019 <- read.csv("data/table2_2019.csv", header = TRUE)
table3.2019 <- read.csv("data/table3_2019.csv", header = TRUE)

table1.2020 <- read.csv("data/table1_2020.csv", header = TRUE)
table2.2020 <- read.csv("data/table2_2020.csv", header = TRUE)
table3.2020 <- read.csv("data/table3_2020.csv", header = TRUE)


### egos


table1.both <- union(table1.2019, table1.2020)
table1.both[is.na(table1.both)] <- 0
table2.both <- union(table2.2019, table2.2020)


# Get egos who stopped using PrEP

egos_stop <- c()
egos_ever<- c()

for (ego in unique(table1.both$egoID) ) {
k_0 <- 0
idx_tmp <- which(table1.both$egoID == ego)
for (i in idx_tmp) {
  prep <- table1.both$PrEP[i]
  k_1 <- k_0
  if(prep == 1 | prep == 2 | prep == 3){k_0 <- 1
  egos_ever <- c(egos_ever, ego)} 
  else{k_0 <- 0}
  if(k_0-k_1 <0){
    egos_stop <- c(egos_stop, ego)
    cat(paste0('ego: ',ego,', visit: ', table1.both$visit[i] ,', '))}
}
}

egos_stop <- unique(egos_stop)
egos_ever <- unique(egos_ever)

table_stop <- table1.both[table1.both$egoID %in% egos_stop,]
alters_stop <- table2.both[table2.both$egoID %in% egos_stop,]

g_ever_stop <- alters_stop %>% group_by(egoID) %>% count() %>% 
  ggplot(aes(x = n)) +
  geom_histogram() + 
  labs(x = 'degree') + 
  labs(y = '# nodes that ever stopped prep')


g_ever_stop
ggsave(g_ever_stop, file = 'Figures/stop_prep.eps', device = 'eps')


### Preparing egos

egos0 <- table1.both[table1.both$visit==1,]
egos0 <- transform(egos0,
                   
                   agecat = as.character(sapply(egos0$age, function(a){
                     if(a<30){"age.gp1"}
                     else if(30<=a & a<40){"age.gp2"}
                     else if(40<=a & a<90){"age.gp3"}
                     #else if(a>=45){"age.gp4"}
                   }))
)

egos0 <- egos0 %>% transform(ethcat = ifelse(white == 1, "white", ifelse(black==1, "black", "other")  ))


egos1 <- table1.both %>% mutate(HIV = ifelse(HIV == 2, 1, 0)) %>% mutate(ART = ifelse(ART ==1 | ART == 2 | ART == 3, 1, 0 )) %>%
  mutate(Group_sex = ifelse(Group_sex>0 & Group_sex<100, 1,0)) %>% mutate(Trade_sex = ifelse(Trade_sex==1,1,0)) %>%
  mutate(PrEP =  ifelse(PrEP == 1 | PrEP == 2 | PrEP ==3, 1, 0 ) )  %>% group_by(egoID) %>% summarise_at(c("HIV","ART","PrEP","Trade_sex", "Group_sex"), sum) %>% 
  mutate(HIV = ifelse(HIV >0, "Positive", "other")) %>% mutate(ART = ifelse(ART > 0, "yes","no")) %>% 
  mutate(PrEP = ifelse(PrEP > 0, "Tr", "Notr")) %>% mutate(Trade_sex = ifelse(Trade_sex >0, "Trade","NoTrade")) %>%
  mutate(Group_sex = ifelse(Group_sex>0, "Group","NoGroup") )


egos <- inner_join(egos0[,c("egoID", "agecat", "ethcat")], egos1, by = "egoID" )

egos  <- transform(egos, treat = egos$PrEP )   
egos <- egos[,c("egoID", "agecat", "ethcat", "treat","Trade_sex", "Group_sex")]
egos <- egos %>% mutate(agecat = as.character(agecat)) %>% mutate(ethcat = as.character(ethcat)) %>% 
  mutate(treat = as.character(treat)) %>% mutate(Trade_sex = as.character(Trade_sex)) %>% mutate(Group_sex = as.character(Group_sex))  

write.csv(egos, "egos_both.csv")



### Preparing alters


table2.both <- union(table2.2019, table2.2020, by = "alterID")
alters10 <- table2.both[-which(is.na(table2.both$age) | (table2.both$age>90) ),] 

alters1 <- alters10[!duplicated(alters10$alterID),]
alters1 <- transform(alters1,
                     agecat = as.character(sapply(alters1$age, function(a){
                       if(a<30){"age.gp1"}
                       else if(30<=a & a<40){"age.gp2"}
                       else if(40<=a & a<90){"age.gp3"}
                       #else if(a>=45){"age.gp4"}
                     }))
)
alters1 <- alters1 %>% transform(ethcat = ifelse(white == 1, "white", ifelse(black==1, "black", "other")  ))
alters1 <- alters1[,c("egoID","alterID","agecat","ethcat")]

table3.both <- union(table3.2019, table3.2020)
table3.both[is.na(table3.both)] <-0
table3.both <- table3.both[which(table3.both$alterID != " "),]

alters2 <- table3.both %>% mutate(HIV = ifelse(HIV == 2, 1, 0)) %>% mutate(ART = ifelse(ART ==1 | ART == 2 | ART ==3, 1, 0 )) %>%
  mutate(Group_sex = ifelse(Group_sex>0 & Group_sex<100, 1,0)) %>% mutate(Trade_sex = ifelse(Trade_sex==1,1,0)) %>%
  mutate(PrEP =  ifelse(PrEP == 1 | PrEP == 2 | PrEP == 3, 1, 0 ) )  %>% group_by(alterID) %>% summarise_at(c("HIV","ART","PrEP","Trade_sex", "Group_sex"), sum) %>% 
  mutate(HIV = ifelse(HIV >0, "Positive", "other")) %>% mutate(ART = ifelse(ART > 0, "yes","no")) %>% 
  mutate(PrEP = ifelse(PrEP > 0, "Tr", "NoTr")) %>% mutate(Trade_sex = ifelse(Trade_sex >0, "Trade","NoTrade")) %>%
  mutate(Group_sex = ifelse(Group_sex>0, "Group","NoGroup") )

alters.trt  <- transform(alters2, treat = alters2$PrEP )   

alters <- inner_join(alters1, alters.trt, by = "alterID")
alters <- alters[,c("egoID", "agecat", "ethcat", "treat", "Trade_sex", "Group_sex")]
alters <- alters %>% mutate(agecat = as.character(agecat)) %>% mutate(ethcat = as.character(ethcat)) %>% 
  mutate(treat = as.character(treat)) %>% mutate(Trade_sex= as.character(Trade_sex)) %>% 
  mutate(Group_sex = as.character(Group_sex))

egos <- egos[which(egos$egoID %in% alters$egoID),]
write.csv(egos, "egos_both.csv")

# Defining sample weights

library(anesrake)

egos <- read.csv("egos_both.csv",  stringsAsFactors=T)

# https://censusreporter.org/profiles/05000US39049-franklin-county-oh/
# "black" "other" "white"
ethtarget <- c(.23, .15, .62)
names(ethtarget) <- c( "black", "other", "white")

agetarget <- c(.19, .16, 0.39)/0.74
names(agetarget) <- c("age.gp1", "age.gp2", "age.gp3")

## compare sample and target distributions

tmp_age <- egos %>% group_by(agecat) %>% count() %>% mutate(percent_sample = n/236) 
tmp_age$percent_tar <- agetarget
tmp_age
#write.csv(tmp_age, "age_weights_both.csv")


tmp_eth <- egos %>% group_by(ethcat) %>% count() %>% mutate(percent_sample = n/236) 
tmp_eth$percent_tar <- ethtarget
tmp_eth
#write.csv(tmp_eth, "eth_weights_both.csv")


targets <- list(ethtarget, agetarget)

names(targets) <- c("ethcat", "agecat")

outsave <- anesrake(targets, egos, caseid=egos$egoID,
                    verbose=TRUE)

caseweights <- data.frame(cases=outsave$caseid, weights=outsave$weightvec)

egos$egoWt <- caseweights$weights

summary(caseweights)
summary(outsave)

egos <- egos %>% transform(agecat = as.character(agecat)) %>% transform(ethcat = as.character(ethcat)) %>%
  transform(treat = as.character(treat)) %>% transform(Trade_sex = as.character(Trade_sex)) %>%
  transform(Group_sex = as.character(Group_sex)) 


NEST <-as.egodata(egos,alters=alters,egoWt = egos$egoWt, egoIDcol="egoID")

write.csv(alters, "alters_both.csv")
write.csv(egos, "egos_both.csv")

save(NEST, file="NESTegodata_both.RData")


# Getting attributes of egos

alters_n <- alters %>% group_by(egoID) %>% count()
alters_age <- alters %>% group_by(egoID, agecat) %>% count() %>% spread(agecat,n) %>% mutate_all(~replace(.,is.na(.),0))
alters_eth <- alters %>% group_by(egoID, ethcat) %>% count() %>% spread(ethcat,n) %>% mutate_all(~replace(.,is.na(.),0))
alters_treat <- alters %>% group_by(egoID, treat) %>% count() %>% spread(treat,n) %>% mutate_all(~replace(.,is.na(.),0))
alters_trade <- alters %>% group_by(egoID, Trade_sex) %>% count() %>% spread(Trade_sex,n) %>% mutate_all(~replace(.,is.na(.),0))
alters_group <- alters %>% group_by(egoID, Group_sex) %>% count() %>% spread(Group_sex,n) %>% mutate_all(~replace(.,is.na(.),0))

alters_des0 <- alters_n %>% left_join(., alters_age, by="egoID") %>% left_join(., alters_eth, by="egoID") %>% left_join(., alters_treat, by="egoID") %>%
  left_join(., alters_trade, by="egoID") %>% left_join(., alters_group, by="egoID")
    
descr_df <- inner_join(egos[c("egoID", "agecat", "ethcat", "treat", "Trade_sex", "Group_sex")], alters_des0, by= "egoID" )

names(descr_df)[7] <- "degree"
dd <- descr_df %>% group_by(degree) %>% count() %>% mutate(Percent = round(n/236,5))
write.csv(descr_df, "sims/attributes/summary_egocentric.csv")
write.csv(dd, "sims/attributes/degree_egocentric.csv")
