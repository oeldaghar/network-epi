# for a given network - find and detect communitites then save cmap for that network

#networkx==3.2.1
import networkx as nx
import os
import csv 

#from more-graph repository
def read_smat(filename, set_undirected = True):
    """ Load an SMAT file into a list-of-lists graph representation. """
    f = open(filename, 'rU')
    hdr = f.readline()
    parts = hdr.split()
    nverts = int(parts[0])
    ncols = int(parts[1])
    nedges = int(parts[2])
    nremedges = nedges # number of remaining edges
    if nverts != ncols:
        raise ValueError(
            'read_smat line 1: requires nrows (%i) = ncols (%i) for graph'%(
            nverts, ncols))
    graph = [ [] for _ in range(nverts) ]
    for lineno, line in enumerate(f):
        parts = line.split()
        if len(parts) == 0: continue
        if nremedges==0: 
            raise ValueError(
                'read_smat line %i: more than %i edges found'%(
                lineno+2, nedges))
        src = int(parts[0])
        dst = int(parts[1])
        if src < 0 or src >= nverts or dst < 0 or dst >= nverts:
            raise ValueError(
                'read_smat line %i: out-of-range edge (%i,%i) found (nverts=%i)'%(
                lineno+2, src, dst, nverts))
        graph[src].append(dst)
        nremedges -= 1
    
    return graph


def get_communities(gpath):
    #load graph
    g = read_smat(gpath)
    g = nx.Graph(dict(enumerate(g)))
    Gcc = sorted(nx.connected_components(g), key=len, reverse=True)
    g = g.subgraph(Gcc[0])

    #get comms 
    comms = nx.community.louvain_communities(g,seed=0)
    cmap = dict() #node_id:comm
    for (i,comm) in enumerate(comms):
        for node in comm:
            cmap[node] = i
    return cmap  



## testing partitions
fpath = "/p/mnt/scratch/network-epi/input/graphs/"
gname = "study-20-draft-150.smat"
gpath = fpath + gname
g = read_smat(gpath)
g = nx.Graph(dict(enumerate(g)))
Gcc = sorted(nx.connected_components(g), key=len, reverse=True)
g = g.subgraph(Gcc[0])


#get comms 
comm_dendogram = nx.community.louvain_partitions(g,seed=0,resolution=0.01)
# for (j,comm_partition) in enumerate(comm_dendogram):
#     print(f"working on level {j}")
#     cmap = dict() #node_id:comm
#     for (i,comm) in enumerate(comm_partition):
#         for node in comm:
#             cmap[node] = i
    
#     cmap = list(cmap.items())
#     cmap.sort(key = lambda x: x[0])
#     cmap = [[x[1]] for x in cmap]


#     fname = os.path.join("/p/mnt/scratch/network-epi/pipeline/data/",gname[:-5],f"cmap-louvain-{j}-0-{gname[:-5]}.csv")
#     with open(fname, 'w') as csvfile:  
#         writer = csv.writer(csvfile)   
#         writer.writerows(cmap)


comms = next(enumerate(comm_dendogram))
[len(x) for x in comms[1]]


dir(igraph)