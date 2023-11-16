#testing epidemic diffusion code
using Random

#allows for execution from command line as well as an ide so files can be run modularlly
mainDir = joinpath("/",split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))]...)

include(joinpath(mainDir,"code","graph-io.jl"))
include(joinpath(mainDir,"code","fast-diffusion.jl"))

#load test graph 
gname = "study-20-draft-150.smat"
gpath = joinpath(mainDir,"input/graphs/")

A = loadGraph(gname,gpath);
nnodes = size(A,1)

#initialize epidemic w/ no quarantining
E = EventEpidemic.SEIRData(A);
seed_node = 227

Random.seed!(11)
l,E = EventEpidemic.epidemic(E,seed_node);

#snodes consistent
total_suscepitble = sum(E.snodes)
total_suscepitble<=nnodes

total_infs = nnodes - total_suscepitble

count(E.itime.==typemax(Int))==total_suscepitble
count(E.ctime.==typemax(Int))==total_suscepitble
count(E.rtime.==typemax(Int))==total_suscepitble

#total_infs instead of nnodes unifected nodes will have times equal to typemax 
sum(E.itime.<E.ctime)==total_infs 

#exposure time consistent with infection time 
sum(E.itime.<=E.ctime)
length(findall(E.itime.>E.ctime))==0

#recovery time consistent with infectious time 
sum(E.ctime.<=E.rtime)
length(findall(E.itime.>E.ctime))==0

#no quarantines 
sum(E.qtimes.==typemax(Int))==nnodes

"""
running the same tests but with quarantining
""" 
E = EventEpidemic.SEIRData(A,qcap=100);
seed_node = 227

Random.seed!(11)
l,E = EventEpidemic.epidemic(E,seed_node,q=true);

#snodes consistent
total_suscepitble = sum(E.snodes)
total_suscepitble<=nnodes

total_infs = nnodes - total_suscepitble

count(E.itime.==typemax(Int))==total_suscepitble
count(E.ctime.==typemax(Int))==total_suscepitble
count(E.rtime.==typemax(Int))==total_suscepitble

#total_infs instead of nnodes unifected nodes will have times equal to typemax 
sum(E.itime.<E.ctime)==total_infs 

#exposure time consistent with infection time 
sum(E.itime.<=E.ctime)
length(findall(E.itime.>E.ctime))==0

#recovery time consistent with infectious time 
sum(E.ctime.<=E.rtime)
length(findall(E.itime.>E.ctime))==0

#quarantines 
#E.qnodes can disagree with E.qtimes in non-zeros. 
#E.qtimes[v] can change when node v is scheduled to be quarantined, even if unsuccessful
#E.qtimes can also be reset if we lack space 

#E.qnodes is true for nodes that were qurantined while infected but reset if Q_S → S
#so if E.qnodes=1, node was quarantine while it was exposed or infected, hence should have been infected 
infected_q_nodes = findall(E.qnodes)

#nodes were infected 
sum(E.snodes[infected_q_nodes].==0)==length(infected_q_nodes)

#consistent with other timings 
sum(E.itime[infected_q_nodes].<typemax(Int))==length(infected_q_nodes)
sum(E.ctime[infected_q_nodes].<typemax(Int))==length(infected_q_nodes)
sum(E.rtime[infected_q_nodes].<typemax(Int))==length(infected_q_nodes)

"""
quarantine reduces epidemic size (in expectation)
when fixing the same random seed and comparing a single diffuion,
random noise can distort this since random numbers are used in a different order
"""

ntrials = 1000
no_quarantine_infs = 0
quarantine_infs = 0
Random.seed!(1)

random_seeds = unique(abs.(rand(Int,ntrials)))
for trial=1:lastindex(random_seeds)
    seed_node = rand(1:nnodes)
    Random.seed!(random_seeds[trial])
    
    #no q 
    update_E!(E,E.beta,E.gamma,0)
    l,E = EventEpidemic.epidemic(E,seed_node,rseed=random_seeds[trial])
    no_quarantine_infs += nnodes-sum(E.snodes)

    #with q 
    update_E!(E,E.beta,E.gamma,100)
    l,E = EventEpidemic.epidemic(E,seed_node,q=true,rseed=random_seeds[trial]);
    quarantine_infs += nnodes-sum(E.snodes)
end

quarantine_infs<no_quarantine_infs

"""
larger beta increases total infections 
"""

ntrials = 10000
beta1,beta2 = 1e-2,1e-2+1e-5
beta1,beta2 = extrema((beta1,beta2))

beta1_infs = 0
beta2_infs = 0

Random.seed!(1)
random_seeds = unique(abs.(rand(Int,ntrials)))
for trial=1:lastindex(random_seeds)
    seed_node = rand(1:nnodes)
    
    #beta 1
    update_E!(E,beta1,E.gamma,0)
    l,E = EventEpidemic.epidemic(E,seed_node,rseed=random_seeds[trial]);
    beta1_infs += nnodes-sum(E.snodes)

    #beta2 
    update_E!(E,beta2,E.gamma,0)
    l,E = EventEpidemic.epidemic(E,seed_node,rseed=random_seeds[trial]);
    beta2_infs += nnodes-sum(E.snodes)
end

beta1_infs<=beta2_infs


"""
smaller gamma increases total infections 
"""

ntrials = 10000
gamma1,gamma2 = 5e-2,5e-2-1e-3
gamma1,gamma2 = extrema((gamma1,gamma2))

gamma1_infs = 0
gamma2_infs = 0

Random.seed!(1)
random_seeds = unique(abs.(rand(Int,ntrials)))
for trial=1:lastindex(random_seeds)
    seed_node = rand(1:nnodes)
    
    #gamma 1
    update_E!(E,E.beta,gamma1,0)
    l,E = EventEpidemic.epidemic(E,seed_node,rseed=random_seeds[trial]);
    gamma1_infs += nnodes-sum(E.snodes)

    #gamma2
    update_E!(E,E.beta,gamma2,0)
    l,E = EventEpidemic.epidemic(E,seed_node,rseed=random_seeds[trial]);
    gamma2_infs += nnodes-sum(E.snodes)
end

#gamma1 smaller so should have more infs 
gamma1_infs>=gamma2_infs
