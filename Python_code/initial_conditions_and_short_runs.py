import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
from collections import Counter
from multiprocessing import Process

import random

#### change if needed

T = 10 # duration of the simulations

nsim = 2
i0_start = 1 # first sampled network
i0_final = 2 # last sampled network


alpha0 = 0.1  # baseline treatment rate
beta0 = 0.05   # transmission rate
scenario0 = '_prep_ic_prep' # the initial condition is PrEP and we consider the PrEP scenario where PrEP nodes have high treatment rate

IC0_start = 1 # first initial condition
IC0_final = 10 # last initial condition

nruns = 10

# With this function, we create a networkx net 

def df_to_graph(E,Deltav):
    Gt = nx.empty_graph()
    N = E.shape[0]
    for i in range(0,N):
        Gt.add_edge(E.iloc[i,0],E.iloc[i,1])
    for i in range(0,len(Deltav)):
        if (i+1) in Gt.nodes():
            Gt.nodes[i+1]['delta'] = Deltav[i]
            Gt.nodes[i+1]['old_label'] = i+1
    G = nx.convert_node_labels_to_integers(Gt,first_label=1)
    
    G = G.to_directed()
    


    N = len(G.nodes)
    A = np.zeros((N,N))

    for edge in G.edges():
        i1 = int(edge[0])
        i2 = int(edge[1])  
        A[i2-1,i1-1] = 1 
        
    deltav = []
    for i in range(1,N+1):
        deltav.append(G.nodes[i]['delta'])
    return G, A, np.array(deltav)


# With the following lines of code, we randomly choose the initial conditions

np.random.seed(seed = 123)
E0_pd = pd.DataFrame()

for i0 in np.arange(1,nsim +1):

    E = pd.read_csv('sims/edges/edges_lcc_sim_'+str(i0)+'.csv')[['V1','V2']]
    deltavM = pd.read_csv('sims/attributes/delta_lcc_sim_'+str(i0)+'.csv')[['V1']]
    summary_df = pd.read_csv('sims/attributes/summary_lcc_sim'+str(i0)+'.csv')
    
    data_comm = pd.read_csv('infomap/attributes/nodes_comm_sim_'+str(i0)+'.csv')

    Deltav = np.array(deltavM.iloc[:,0])

    G0,A0, deltav0 = df_to_graph(E,Deltav) # careful: the labels of the largest connected component are changed here!
    
    gamma = 1

    ICs =[]
    for node in G0.nodes():
        ICs_large = []
        ICs_small = []
        for node in G0.nodes():
            if (data_comm.loc[node-1,'type_comm_tr'] == 'large') and (data_comm.loc[node-1, 'type_comm_no_tr'] == 'large') and (np.round(G0.nodes[node]['delta'],2) == gamma):
                ICs_large.append(node)
            elif (data_comm.loc[node-1,'type_comm_tr'] == 'small') and (data_comm.loc[node-1, 'type_comm_no_tr'] == 'large') and (np.round(G0.nodes[node]['delta'],2) == gamma):
                ICs_small.append(node)

    ICs = ICs + list(np.random.choice(ICs_large,size = 50, replace = False))
    if len(ICs_small) < 50:
        ICs = ICs + list(np.random.choice(ICs_small,size = 50, replace = True))
    else:
        ICs = ICs + list(np.random.choice(ICs_small,size = 50, replace = False))
            
    E0_pd['node'] = list(np.array(ICs).reshape(1,-1)[0])
    E0_pd['old_label'] = [G0.nodes[node]['old_label'] for node in list(np.array(ICs).reshape(1,-1)[0])]
    E0_pd['delta'] = [G0.nodes[node]['delta'] for node in list(np.array(ICs).reshape(1,-1)[0])]
    E0_pd['degree'] =[len(list(G0.neighbors(node))) for node in list(np.array(ICs).reshape(1,-1)[0])]
    E0_pd['PrEP_comm'] = list(data_comm.loc[list(np.array(ICs).reshape(1,-1)[0]-1),'type_comm_tr'])
    E0_pd.to_csv('dynamics/initial_conditions/initial_conditions_sim_' + str(i0) +'.csv')



### DISEASE DYNAMICS: With this function, we execute a single run of our syphilis model

def syphsimT(G, TargetN, E0, beta, alpha):
    
    
    alpha1 = alpha + pabx + 1/1.5*(pt1) # I1 ->T1
    alpha2 = alpha + pabx + 1/3.6*(pt2) # I2 -> T1
    alpha3 = alpha + pabx + 1/6.9*(pt3) # L1 -> T2
    alpha4 = alpha + pabx # L2 -> T3


    alpha1p = alpha + pabx + 1/1.5*(pt1) + 1/1.5 # I1 ->T1 in PrEP
    alpha2p = alpha + pabx + 1/3.6*(pt2) + 1/1.5 # I2 -> T1 in PrEP
    alpha3p = alpha + pabx + 1/6.9*(pt3) + 1/1.5 # L1 -> T2 in PrEP
    alpha4p = alpha + pabx + 1/1.5 # L2 -> T3 in PrEP
    
    Nodes = list(G.nodes())

    G.nodes[E0]['state'] = 'E'
    
    t = 0
    fs = 1
    tv = []
    Nv = []
    node_i = E0
    
    while t < T:
        
        
        SIl = []
        Sl = []
        El = []
        I1l = []
        I2l = []
        T1l = []
        L1l = []
        L2l = []
        T2l = []
        T3l = []
        I1pl = []
        I2pl = []
        L1pl = []
        L2pl = []
        
        for node in list(G.nodes()):
            if G.nodes[node]['state'] == 'S':
                Sl.append(node)
            if G.nodes[node]['state'] == 'E':
                El.append(node)
            if G.nodes[node]['state'] == 'I1':
                I1l.append(node)
            if G.nodes[node]['state'] == 'I2':
                I2l.append(node)
            if G.nodes[node]['state'] == 'T1':
                T1l.append(node)
            if G.nodes[node]['state'] == 'L1':
                L1l.append(node)
            if G.nodes[node]['state'] == 'L2':
                L2l.append(node)
            if G.nodes[node]['state'] == 'T2':
                T2l.append(node)
            if G.nodes[node]['state'] == 'T3':
                T3l.append(node)
            if G.nodes[node]['state'] == 'I1p':
                I1pl.append(node)
            if G.nodes[node]['state'] == 'I2p':
                I2pl.append(node)
            if G.nodes[node]['state'] == 'L1p':
                L1pl.append(node)
            if G.nodes[node]['state'] == 'L2p':
                L2pl.append(node)
                
        for node in Sl:
            for nodep in G.neighbors(node):
                nst = G.nodes[nodep]['state']
                if nst == 'I1' or nst == 'I2' or nst == 'I1p' or nst == 'I2p':
                    SIl.append([node,nodep])
                    
        nSI = len(SIl); nS = len(Sl); nE = len(El); nI1 = len(I1l); nI2 = len(I2l)
        nT1 = len(T1l); nT2 = len(T2l); nT3 = len(T3l); nL1 = len(L1l); nL2 = len(L2l)
        nI1p = len(I1pl); nI2p = len(I2pl); nL1p = len(L1pl); nL2p = len(L2pl)
        
        
        tv.append(t)
        Nv.append([t,nS,nE,nI1,nI1p,nI2,nI2p,nL1,nL1p,nL2,nL2p,nT1,nT2,nT3, node_i])
        
        lamSE = beta*nSI; lamEI = delt* nE; lamI1T1 = alpha1*nI1; lamI1I2 = gamma1*nI1
        lamI2T1 = alpha2*nI2; lamI2L1 = gamma2*nI2; lamT1S = lam1*nT1; lamL1T2 = alpha3*nL1
        lamL1L2 = gamma3*nL1; lamL2T3 = alpha4*nL2; lamT2S = lam2*nT2; lamT3S = lam3*nT3
        lamI1pI2p = gamma1*nI1p ; lamI2pL1p = gamma2*nI2p; lamL1pL2p = gamma3*nL1p 
        lamI1pT1 = alpha1p*nI1p; lamI2pT1 = alpha2p*nI2p; lamL1pT2 = alpha3p*nL1p; lamL2pT3 = alpha4p*nL2p
        
        lamv = np.array([lamI1T1, lamI1I2, lamI2T1, lamI2L1, lamT1S, lamL1T2, lamL1L2, lamL2T3, lamT2S, lamT3S, lamI1pI2p, lamI2pL1p, lamL1pL2p, lamI1pT1, lamI2pT1, lamL1pT2, lamL2pT3])
        
        if nE ==0 and nI1 == 0 and nI2==0 and nI1p==0 and nI2p==0:
            break
        statepop = ['T1','I2','T1','L1','S','T2','L2','T3','S','S','I2p','L1p','L2p','T1','T1','T2','T2']
        listpop =  [I1l, I1l, I2l, I2l, T1l, L1l, L1l, L2l, T2l, T3l,I1pl,I2pl,L1pl, I1pl,I2pl, L1pl,L2pl]
        
        lamTot = np.sum(lamv) + lamSE + lamEI
        t = t + np.random.exponential((1/lamTot),1)[0]
        u = random.uniform(0,1)
        
        
        if u < lamSE/lamTot:
            pair = random.choice(SIl)
            nodeS = pair[0]
            node_i = nodeS
            G.nodes[nodeS]['state'] = 'E'
            fs = fs+1
        
        if lamSE/lamTot < u and u < (lamSE+lamEI)/lamTot:
            nodeInf = random.choice(El)
            node_i = nodeInf
            if nodeInf in TargetN:
                G.nodes[nodeInf]['state'] = 'I1p'
            else:
                G.nodes[nodeInf]['state'] = 'I1'
        
        for i in range(0,len(statepop)):
            if (lamSE + lamEI + np.sum(lamv[:i]))/lamTot < u and u < (lamSE+lamEI+np.sum(lamv[:(i+1)]))/lamTot:
                nodeCh = random.choice(listpop[i])
                node_i = nodeCh
                G.nodes[nodeCh]['state'] = statepop[i]
    
    for node in list(G.nodes()):
        G.nodes[node]['state'] = 'S'
    return tv, Nv
        

# Parameters of the syphilis model

# Proportion seeking treatment due to infection 
pt1 = 0.15
pt2 = 0.25
pt3 = 0.25

pabx = 0.001/12

delt = 1/0.9

lam1 = 1/0.25 # T1 -> S
lam2 = 1/0.25 # T2 -> S
lam3 = 1/60 # T3 -> S

gamma1 = 1/1.5*(1-pt1) # 1/1.5 ?  I1 -> I2
gamma2 = 1/3.6*(1-pt2) # 1/3.6 ? I2 -> L1
gamma3 = 1/6.9*(1-pt3) # 1/6.9 ? L1 -> L2


# With the following function we run our syphilis model on sample networks of our ERGM        

def runsim1(beta, i_start, i_final, k1, k2, scenario, alpha):

    for i0 in range(i_start, i_final + 1):
    
        E = pd.read_csv('sims/edges/edges_lcc_sim_'+str(i0)+'.csv')[['V1','V2']]
        deltavM = pd.read_csv('sims/attributes/delta_lcc_sim_'+str(i0)+'.csv')[['V1']]
        summary_df = pd.read_csv('sims/attributes/summary_lcc_sim'+str(i0)+'.csv')
        Deltav = np.array(deltavM.iloc[:,0])
        E0s = list(pd.read_csv('dynamics/initial_conditions/initial_conditions_sim_' + str(i0) +'.csv')['node'])

        G0,A0, deltav0 = df_to_graph(E,Deltav)


        Non_treat = []
        Treat = [] # this is gonna be TargetN

        for node in G0.nodes():
            if G0.nodes[node]['delta'] != 1:
                Non_treat.append(node)    

        for node in G0.nodes():
            if G0.nodes[node]['delta'] == 1:
                Treat.append(node)

    
        if scenario == '_prep_ic_prep':
            TargetN = Treat
            #IC = Treat

        if scenario == '_prep_ic_non_prep':
            TargetN = Treat
            #IC = Non_treat

        if scenario == '_non_prep':
            TargetN = []
            #IC = list(G0.nodes())

        for node in list(G0.nodes()):
            G0.nodes[node]['state'] = 'S'

            
        for k in range(k1,k2+1):
            E0 = E0s[k-1]
            for kk in range(1,nruns+1):
                tv, Nv = syphsimT(G0, TargetN, E0 , beta, alpha)
                Nv_df = pd.DataFrame(Nv, columns = ['t','nS','nE','nI1','nI1p','nI2','nI2p','nL1','nL1p','nL2','nL2p','nT1','nT2','nT3', 'node'])
                Nv_df.to_csv('dynamics/short_runs/short_IC_'+str(k)+'_run_'+str(kk)+'_sim_'+str(i0)+'.csv')


runsim1(beta0,i0_start,i0_final, IC0_start, IC0_final, scenario0, alpha0)

# Take attributes from short runs

net_long = []
scenario_long =[]
duration = []
small_large = []
equilibria = []
degree= []
year = []
node_long = []
delta_long = []
PrEP_comm_long = []
no_PrEP_comm_long = []

for i0 in range(1,nsim+1):
    
    E = pd.read_csv('sims/edges/edges_lcc_sim_'+str(i0)+'.csv')[['V1','V2']]
    deltavM = pd.read_csv('sims/attributes/delta_lcc_sim_'+str(i0)+'.csv')[['V1']]
    nodes_comm = pd.read_csv('infomap/attributes/nodes_comm_sim_'+str(i0)+'.csv')
    Deltav = np.array(deltavM.iloc[:,0])

    G0,A0, deltav0 = df_to_graph(E,Deltav) 
    
    for scenario in ['_prep_ic_prep']:
        for IC in range(IC0_start,IC0_start+1):
            for k in range(1,nruns+1):
                Nv = pd.read_csv('dynamics/short_runs/short_IC_'+str(IC)+'_run_'+str(k)+'_sim_'+str(i0)+'.csv') #[t,nS,nE,nI1,nI1p,nI2,nI2p,nL1,nL1p,nL2,nL2p,nT1,nT2,nT3]
                E0 = Nv.node.iloc[0]
                
                scenario_long.append(scenario)
                duration.append(Nv.t.iloc[-1])
                degree.append(len(list(G0.neighbors(E0))))
                delta_long.append(np.round(G0.nodes[E0]['delta'],2))
                node_long.append(E0)
                PrEP_comm_long.append(nodes_comm.type_comm_tr.iloc[E0-1])
                no_PrEP_comm_long.append(nodes_comm.type_comm_no_tr.iloc[E0-1])
                net_long.append(i0)
                
                Ny = (Nv.nI1+Nv.nI2+Nv.nI1p+Nv.nI2p)*100/len(G0.nodes())

                fs_tmp = list(Ny[np.array(Nv.t)>T-40])
                if fs_tmp == []:
                    equilibria.append(0)
                else: 
                    equilibria.append(np.mean(fs_tmp))

                if Nv.t.iloc[-1] < 6:
                    small_large.append('small')

                if Nv.t.iloc[-1] >6:
                    small_large.append('large')


output_pd = pd.DataFrame({'network':net_long, 'scenario': scenario_long, 'duration':duration,'outbreak_type':small_large,'equilibria':equilibria, 'E0': node_long, 'delta':delta_long, 'degree':degree, 'PrEP_comm_type':PrEP_comm_long, 'no_PrEP_comm_type':no_PrEP_comm_long})
output_pd.to_csv('dynamics/summary_all_short_runs.csv')
