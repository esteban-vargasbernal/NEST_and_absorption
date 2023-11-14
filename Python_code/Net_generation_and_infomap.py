import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from numpy.linalg import inv
from scipy.linalg import expm, sinm, cosm
import infomap


### change these if needed

n_sim = 2 # number of sampled networks
t_0 = 0.3 # Markov time for which we take a partition


# With this function, we get a networkx net with absorption as a node attribute


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


### We extract the largest connected component of each sampled network from our ERGM

deltavM = pd.read_csv('../R_code/sims/attributes/delta_vector.csv')[['V1']]
Deltav = np.array(deltavM.iloc[:,0])

# Defining the two possible node-absorption rates

gamma_1 = 1/154/(1/154 + 1/45) # 46 days in I1, 108 days in I2, 45 days if PrEP 
gamma_2= 1
Deltav[Deltav==0.2] = gamma_1
Deltav[Deltav==1 ] = gamma_2

largest_size_cc_list = []
largest_cc_list = []

for i in range(1,n_sim+1):

    E = pd.read_csv('../R_code/sims/edges/edges_sim_'+str(i)+'.csv')[['V1','V2']]
    summary = pd.read_csv("../R_code/sims/attributes/summary_sim"+str(i)+".csv")
    G0_start,A0, deltav0 = df_to_graph(E,Deltav) # this graph does not have nodes with degree zero

    cc_ls = list(nx.connected_components(G0_start.to_undirected()))
    size_cc = []

    for cc in cc_ls:
        size_cc.append(len(cc))

    largest_size_cc_list.append(np.max(np.array(size_cc)))
    max_idx = np.where(np.array(size_cc) == np.max(np.array(size_cc)))[0][0]
    
    H = G0_start.subgraph(list(cc_ls[max_idx]))
    H = nx.convert_node_labels_to_integers(H,first_label=1)
    H = H.to_directed()
    E_new = pd.DataFrame(np.array(list(H.edges())), columns = ['V1','V2'])
    
    old_labels=[]
    for node in H.nodes():
        old_labels.append(H.nodes[node]['old_label']-1)

    deltav_new = pd.DataFrame(np.array([H.nodes[node]['delta'] for node in H.nodes()]), columns=['V1'])
    summary_new = summary.iloc[old_labels,:]    
    deg_new = pd.DataFrame(summary_new[['node','degree']].sort_values(by = ['degree']))
    
    
    E_new.to_csv('sims/edges/edges_lcc_sim_'+str(i)+'.csv')
    deltav_new.to_csv('sims/attributes/delta_lcc_sim_'+str(i)+'.csv')
    deg_new.to_csv("sims/attributes/degree_lcc_sim_"+str(i)+".csv")
    summary_new.to_csv("sims/attributes/summary_lcc_sim"+str(i)+".csv")
    
large_size = pd.DataFrame(np.array(largest_size_cc_list), columns = ['size'] )
large_size.to_csv('sims/attributes/sizes_connected_comp.csv')


### With the following five functions, we can run InfoMap for absorbing random walks

# With this function, we get a graph Laplacian

def infomap_abs_L(G,A,deltav, typeAdelta, h):
            
    N = A.shape[0]    
    outfl = A.sum(axis=0)
    W = np.diag(outfl)
    Ddelta = np.diag(deltav)
    
    Ltilde = (W-A).dot(inv(h*W + Ddelta)) # use 'Linear', 'Exp'
    LtildeNor = (np.identity(N)-A.dot(inv(W))).dot(inv(Ddelta)) 
           
    if (typeAdelta == 'Linear'): 
        Ltilde = (W-A).dot(inv(h*W + Ddelta))
        tmax = np.multiply((deltav + h*outfl),1/outfl).min()
        
    if (typeAdelta == 'Exp'): 
        Ltilde = (W-A).dot(inv(h*W + Ddelta))
        tmax = 10
        
    if (typeAdelta == 'LinearNor'): 
        Ltilde = (np.identity(N)-A.dot(inv(W))).dot(inv(Ddelta)) 
        tmax = deltav.min()
        
    if (typeAdelta =='ExpNor'): 
        Ltilde = (np.identity(N)-A.dot(inv(W))).dot(inv(Ddelta)) 
        tmax = 10
        
    if typeAdelta == 'NormalInfo':
        Ltilde = A
        tmax = 1

    return Ltilde, tmax

# With this function, we get an adjacency matrix

def infomap_abs_input(Ltilde, t, typeAdelta, h):
    
    N = Ltilde.shape[0]

    if typeAdelta == 'Linear': 
        Adelta = np.identity(N) - t*Ltilde
            
    if typeAdelta == 'Exp':
        Adelta = expm(-t*Ltilde)
        
    if typeAdelta == 'LinearNor': 
        Adelta = np.identity(N) - t*Ltilde
            
    if typeAdelta == 'ExpNor':
        Adelta = expm(-t*Ltilde)
    
    
    if typeAdelta == 'NormalInfo':
        Adelta = Ltilde

    return Adelta

# With this function, we run InfoMap for absorbing random walk for a specific Markov time t

def infomap_abs(G, A, Ltilde,  deltav,  t, typeAdelta, h, Get_Mod):
    
    Adelta = infomap_abs_input(Ltilde, t, typeAdelta,h)

    N = Adelta.shape[0]
    im = infomap.Infomap("--directed --two-level --teleportation-probability 0")


    for i in range(0,N):
        for j in range(0,N):
            im.add_link(i,j,Adelta[j,i])
            im.add_node(i,str(i+1))
            im.add_node(j,str(j+1))
    im.run()

    if Get_Mod == True:
        Modv = pd.DataFrame(columns = ["node","Module"])
        ii = 0
        for node, module in im.modules:
            Modv.loc[ii,:] = np.array([node, module])
            ii = ii+1


        Modv = pd.DataFrame(Modv, columns = ["node","Module"])
        Modv = Modv.sort_values("node")
        cluster = list(Modv.loc[:,"Module"])

        return im.numTopModules(), cluster 

    else:
        return im.numTopModules()

# With this function, we run InfoMap for absorbing random walks with Nt equally spaced Markov times 
# in the interval [t0,tf].  

def InfoAbs(G, A, deltav,  t0, tf, Nt, typeAdelta, h): 
    if typeAdelta != 'NormalInfo':
        
        Ltilde, tmax = infomap_abs_L(G,A,deltav, typeAdelta, h)
        
        Nmodules = []

        if typeAdelta == 'Linear' or typeAdelta == 'LinearNor':
            print('Max. possible value of t: ', tmax)
            
        tend = np.array([tf,tmax]).min()
        th = (tend-t0)/Nt
        tv =np.arange(t0,tend + th,th)

        for t in tv:
            n_tmp = infomap_abs(G, A, Ltilde, deltav, t, typeAdelta, h, False)
            Nmodules.append(n_tmp)


                                           
        return tv, Nmodules
    else:
        print('Number of communities inferred by regular InfoMap: ', infomap_abs(G, deltav, 1, typeAdelta, h)[0])
        

# This function extracts the partitions yielded by InfoMap for absorbing random walks
# for the Markov times in the vector tstar_ls

def InfoAbs_tstar(G, A,  deltav, typeAdelta, h, tstar_ls): 
                
    cluster_ls = []

    for t in tstar_ls:
    
        Ltilde = infomap_abs_L(G,A,deltav, typeAdelta, h)[0]
        _ , cluster = infomap_abs(G,A, Ltilde, deltav, t, typeAdelta,h, True)
        cluster_ls.append(cluster)
                                      
    return  cluster_ls

    

### Generate the partitions for t_0

for i0 in range(1,n_sim + 1):

    E = pd.read_csv('sims/edges/edges_lcc_sim_'+str(i0)+'.csv')[['V1','V2']]
    deltavM = pd.read_csv('sims/attributes/delta_lcc_sim_'+str(i0)+'.csv')[['V1']]
    summary_df = pd.read_csv('sims/attributes/summary_lcc_sim'+str(i0)+'.csv')
    Deltav = np.array(deltavM.iloc[:,0])


    G0,A0, deltav0 = df_to_graph(E,Deltav) 

    tv_N_PrEP = InfoAbs(G0, A0, deltav0, 0.01, 1, 2, 'Exp', 0)

    tv_N_PrEP_df = pd.DataFrame(np.array(tv_N_PrEP).transpose(), columns = ['t','nc'])    
    tv_N_PrEP_df.to_csv('infomap/PrEP/nc_vs_t_PrEP_sim_'+ str(i0) + '.csv')

    cluster_ls_trt = InfoAbs_tstar(G0, A0,  deltav0, 'Exp', 0, [t_0])

    np.savetxt('infomap/PrEP/partition_PrEP_sim_'+str(i0)+'.csv',cluster_ls_trt[0],delimiter = ',')
    
    tv_N_no_PrEP = InfoAbs(G0, A0, np.array([np.min(deltav0)]*len(deltav0)), 0.01, 1, 2, 'Exp', 0)

    tv_N_no_PrEP_df = pd.DataFrame(np.array(tv_N_no_PrEP).transpose(), columns = ['t','nc'])    
    tv_N_no_PrEP_df.to_csv('infomap/no_PrEP/nc_vs_t_no_PrEP_sim_'+ str(i0) + '.csv')

    cluster_ls_non_trt = InfoAbs_tstar(G0, A0,  np.array([np.min(deltav0)]*len(deltav0)) , 'Exp', 0, [t_0])

    
    np.savetxt('infomap/no_PrEP/partition_no_PrEP_sim_'+str(i0)+'.csv',cluster_ls_non_trt[0],delimiter = ',')
    
# Plot number of communities versus Markov time

for i0 in np.arange(1,n_sim + 1):
    t0=1
    tv_N_PrEP = pd.read_csv('infomap/PrEP/nc_vs_t_PrEP_sim_'+ str(i0) + '.csv')
    tv_N_no_PrEP = pd.read_csv('infomap/no_PrEP/nc_vs_t_no_PrEP_sim_'+ str(i0) + '.csv')
    
    tv_PrEP = tv_N_PrEP.t.iloc[t0:]
    nc_PrEP = tv_N_PrEP.nc.iloc[t0:] 
    
    tv_no_PrEP = tv_N_no_PrEP.t.iloc[t0:]
    nc_no_PrEP = tv_N_no_PrEP.nc.iloc[t0:] 
    
    

    plt.plot(tv_PrEP,nc_PrEP, 'r', color = '0.4')
    plt.plot(tv_no_PrEP,nc_no_PrEP, 'g', color = '0.7')
plt.axvline(x=0.6, color = 'black', linestyle = 'dashed')

plt.xlabel(r'$t$')
plt.ylabel('Number of communities')
plt.legend(['Defferent absorptions', 'Same absorptions', 't=0.3'])
plt.savefig('dynamics/Figures/number_of_comm_vs_t_t0_0_6.eps', format = 'eps')

# Reading community attributes and PrEP neighbors
    
for i0 in np.arange(1,n_sim + 1):
    E = pd.read_csv('sims/edges/edges_lcc_sim_'+str(i0)+'.csv')[['V1','V2']]
    deltavM = pd.read_csv('sims/attributes/delta_lcc_sim_'+str(i0)+'.csv')[['V1']]
    cluster_PrEP = np.loadtxt('infomap/PrEP/partition_PrEP_sim_'+str(i0)+'.csv', delimiter = ',')
    cluster_no_PrEP = np.loadtxt('infomap/no_PrEP/partition_no_PrEP_sim_'+str(i0)+'.csv',delimiter = ',')
    
    Deltav = np.array(deltavM.iloc[:,0])

    G0,A0, deltav0 = df_to_graph(E,Deltav) # here there is a relabeling that is consitent with the relabeling for the clusterings
    
    nc1 = len(np.unique(cluster_PrEP))
    sc1 = [len(cluster_PrEP[cluster_PrEP==i]) for i in range(1,nc1+1)]
    count_sc1 = np.array([(j,sc1.count(j)) for j in np.unique(sc1)])
    
    nc2 = len(np.unique(cluster_no_PrEP))
    sc2 = [len(cluster_no_PrEP[cluster_no_PrEP==i]) for i in range(1,nc2+1)]
    count_sc2 = np.array([(j,sc2.count(j)) for j in np.unique(sc2)])
    
    size_max_PrEP = np.max(count_sc1[:,0])
    size_max_no_PrEP = np.max(count_sc2[:,0])
    
    for i in range(1, len(G0.nodes())+1):
        ct = cluster_PrEP
        cnt = cluster_no_PrEP
        G0.nodes[i]['comm_PrEP'] = ct[i-1]
        G0.nodes[i]['comm_no_PrEP'] = cnt[i-1]
        G0.nodes[i]['size_PrEP_comm'] = len(ct[ct == ct[i-1]])
        G0.nodes[i]['size_no_PrEP_comm'] = len(cnt[cnt == cnt[i-1]])

        
    # This section creates a dataframe with the proportion of PrEP nodes inside the PrEP communities that intersect the large no-PrEP community
    
    data_cs = np.array([(G0.nodes[node]['comm_PrEP'], np.round(G0.nodes[node]['delta'],2)) for node in G0.nodes() if (((G0.nodes[node]['size_no_PrEP_comm'] == int(size_max_no_PrEP)) ))])
    data_cs_df = pd.DataFrame(data_cs,columns = ['comm','delta']).sort_values('comm')
    data_cs_df.groupby('comm').count()

    nodes_tmp = []

    for comm in np.unique(data_cs_df.comm):    

        ntr = 0
        size = 0

        for node in G0.nodes():

            if G0.nodes[node]['comm_PrEP'] == comm:
                size = size +1
                if G0.nodes[node]['delta'] == 1:
                    ntr = ntr+1

        nodes_tmp.append((comm, size, ntr, np.round(ntr/size,3)))

    data_comm_small = pd.DataFrame(np.array(nodes_tmp), columns = ['comm','size','n_trt','prop'])        
    data_comm_small.to_csv('infomap/attributes/PrEP_prop_per_comm_inside_large_no_PrEP_sim'+str(i0)+'.csv')
    
    # This section creates a dataframe with the community attributes of each node

    tmp = []
    for node in G0.nodes():
        if G0.nodes[node]['delta'] == 1:
            tmp.append((node,'Tr',len(list(G0.neighbors(node))), G0.nodes[node]['size_PrEP_comm'], G0.nodes[node]['size_no_PrEP_comm'], G0.nodes[node]['old_label']))
        else: 
            tmp.append((node,'NoTr',len(list(G0.neighbors(node))), G0.nodes[node]['size_PrEP_comm'], G0.nodes[node]['size_no_PrEP_comm'], G0.nodes[node]['old_label']))
    data_comm = pd.DataFrame(np.array(tmp), columns = ['node', 'treat', 'degree', 'size_comm_tr', 'size_comm_no_tr','old_label'])

    data_comm = data_comm.assign(type_comm_tr = lambda x: ['large' if (y == str(size_max_PrEP)) else 'small' for y in x.size_comm_tr]).assign(type_comm_no_tr = lambda x: ['large' if (y == str(size_max_no_PrEP)) else 'small' for y in x.size_comm_no_tr])
    data_comm.to_csv('infomap/attributes/nodes_comm_sim_'+str(i0)+'.csv')
    
    
    # This section is for the PrEP neighbors in neighborhoods of radius 4 of PrEP nodes in the large no - PrEP community 
    
    nn =4
    n_trt = []
    n_dg = []
    list_nodes = []
    for node in G0.nodes():
        if G0.nodes[node]['size_no_PrEP_comm'] == int(size_max_no_PrEP) and G0.nodes[node]['size_PrEP_comm'] < int(size_max_PrEP) and G0.nodes[node]['delta']==1:
            list_nodes.append(G0.nodes[node]['old_label'])
            k = 0
            n_dg.append(len(list(G0.neighbors(node))))
            X_tmp = list(nx.ego_graph(G0,node, radius = nn).nodes())
            for node2 in X_tmp:
                if G0.nodes[node2]['delta'] == 1:
                    k = k+1
            n_trt.append(k/len(X_tmp))

    ngb_r4_small = pd.DataFrame({'node': list_nodes,'prop_trt':n_trt, 'degree': n_dg})

    n_trt2 = []
    n_dg2 = []
    list_nodes = []
    for node in G0.nodes():
        if G0.nodes[node]['size_no_PrEP_comm'] == int(size_max_no_PrEP) and G0.nodes[node]['size_PrEP_comm'] == int(size_max_PrEP) and G0.nodes[node]['delta']==1:
            list_nodes.append(G0.nodes[node]['old_label'])
            k = 0
            n_dg2.append(len(list(G0.neighbors(node))))
            X_tmp = list(nx.ego_graph(G0,node, radius = nn).nodes())
            for node2 in X_tmp:
                if G0.nodes[node2]['delta'] == 1:
                    k = k+1
            n_trt2.append(k/len(X_tmp))

    ngb_r4_large = pd.DataFrame({'node':list_nodes,'prop_trt':n_trt2, 'degree': n_dg2})

    ngb_r4_large.to_csv('infomap/attributes/PrEP_neighbors_large_comm_sim_'+str(i0)+'.csv')
    ngb_r4_small.to_csv('infomap/attributes/PrEP_neighbors_small_comm_sim_'+str(i0)+'.csv')