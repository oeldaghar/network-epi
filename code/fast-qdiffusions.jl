##
using Distributed,ProgressMeter
# addprocs(2)
@everywhere mainDir ="/p/mnt/scratch/network-epi/" 
@everywhere include(joinpath(mainDir,"code/fast-diffusion.jl"))
@everywhere include(joinpath(mainDir,"code/graph-io.jl"))
@everywhere include(joinpath(mainDir,"code/data-io.jl"))

gpath = joinpath(mainDir,"pipeline/graphs/")
parentDst = joinpath(mainDir,"pipeline/data/")

# gs = ["mexico","anony","dblp-cc",
#     "study-11-45","study-11-50",
#     "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5",
#     "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.9-0.0-100-5",
#     "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17.smat",
#     "enron","uf","commutes","astro","penn","ucf",
#     "usf","modberkeley","modfsu","modharvard","modnotre-dame","modstanford",
#     "moduillinois","modunc","modwisconsin"]#,"flickr","livejournal"]

function read_inf_data(gname;dloc=dloc,beta=0.1,gamma=0.05,method="seir",dtype="cinfs")
    if endswith(gname,".smat")
        g = gname[1:end-5]
    else
        g = gname
    end
    fname = dloc*"$dtype-$g-$method-$beta-$gamma.txt"
    if lowercase(dtype)=="cinfs" || lowercase(dtype)=="tinfs"
        res = Vector{Vector{Int}}()
        open(fname,"r") do io
            while !eof(io)
                push!(res,parse.(Int,split(readline(io),",")))
            end
        end
    else
        error("no such file or wrong data format. this reads Vector{Vector{Int}} data")
    end
    return res
end

# bs = [collect(0.001:0.001:0.01);collect(0.02:0.01:0.1)]
qpercents = collect(0:15)
nnodes,ktrials = 50,1
tmax = 30000

# function _checkgraph_old(gname::String,bs::Vector{Float64};dst::String="",gamma::Float64=0.05,method::String="sir",
#                         qpercents::Vector{T}=collect(0:15),nnodes::Int=50) where T<:Union{Int,Float64}
#     """checks to see if we have already run diffusions for this graph"""
#     if endswith(gname,".smat")
#         g = gname[1:end-5]
#     else
#         g = gname
#     end

#     fnames = filter(c->endswith(c,".txt") && occursin(g,c),readdir(dst))
#     betabool = ones(Bool,length(bs)) #indicates if we need to run bs[i] or not
#     for (i,b) in enumerate(bs)
#         if "tinfs-$g-$method-$b-$gamma.txt" in fnames && "cinfs-$g-$method-$b-$gamma.txt" in fnames && "scounts-$g-$method-$b-$gamma.txt" in fnames
#             tinfs = read_inf_data(g,dloc=dst,beta=b,gamma=gamma,method=method,dtype="tinfs")
#             if length(tinfs)==length(qpercents)*nnodes 
#                 betabool[i] = false #files exist and are of right length
#             end
#         end
#     end
#     return betabool
# end

function _checkgraph(gname::String,bs::Vector{Float64};dst::String="",gamma::Float64=0.05,method::String="sir",
                        qpercents::Vector{T}=collect(0:15),nnodes::Int=50) where T<:Union{Int,Float64}
    """checks to see if we have already run diffusions for this graph"""
    if endswith(gname,".smat")
        g = gname[1:end-5]
    else
        g = gname
    end

    fnames = filter(c->endswith(c,".txt") && occursin(g,c),readdir(dst))

    #quick patch to fix issues with betas
    oldbetas = get_betas(gname)
    curr_betas = sort(unique(vcat(oldbetas,bs)))
    #end patch 
        
    betabool = ones(Bool,length(curr_betas)) #indicates if we need to run bs[i] or not
    for (i,b) in enumerate(curr_betas)
        if "tinfs-$g-$method-$b-$gamma.txt" in fnames && "cinfs-$g-$method-$b-$gamma.txt" in fnames && "scounts-$g-$method-$b-$gamma.txt" in fnames
            tinfs = read_inf_data(g,dloc=dst,beta=b,gamma=gamma,method=method,dtype="tinfs")
            if length(tinfs)==length(qpercents)*nnodes 
                betabool[i] = false #files exist and are of right length
            end
        end
    end
    return curr_betas,betabool
end


function _checkgraph_new(gname::String,bs::Vector{Float64};dst::String="",gamma::Float64=0.05,method::String="sir",
                        qpercents::Vector{T}=collect(0:15),nnodes::Int=50,
                        oldbetas::Vector{Float64}=Vector{Float64}()) where T<:Union{Int,Float64}
    """checks to see if we have already run diffusions for this graph"""
    if endswith(gname,".smat")
        g = gname[1:end-5]
    else
        g = gname
    end

    fnames = filter(c->endswith(c,".txt") && occursin(g,c),readdir(dst))

    #quick patch to fix issues with betas
    # oldbetas = get_betas(gname)
    curr_betas = sort(unique(vcat(oldbetas,bs)))
    #end patch 
        
    betabool = ones(Bool,length(curr_betas)) #indicates if we need to run bs[i] or not
    for (i,b) in enumerate(curr_betas)
        if "tinfs-$g-$method-$b-$gamma.txt" in fnames && "cinfs-$g-$method-$b-$gamma.txt" in fnames && "scounts-$g-$method-$b-$gamma.txt" in fnames
            tinfs = read_inf_data(g,dloc=dst,beta=b,gamma=gamma,method=method,dtype="tinfs")
            if length(tinfs)==length(qpercents)*nnodes 
                betabool[i] = false #files exist and are of right length
            end
        end
    end
    return curr_betas,betabool
end


bs = unique(vcat(1e-3:1e-3:1e-2,2e-2:1e-2:1e-1,1.1e-1:1e-2:2e-1))
qpercents = collect(0:15)

# gs = [ "study-11-2023-0-noweak.smat",
#     "study-11-2022-1.smat",
#     # "study-11-2022-10.smat",
#     # "study-11-2022-20.smat",
#     # "study-11-2022-30.smat",
#     # "study-11-2022-40.smat",
#     # "study-11-2022-45.smat",
#     "study-11-2022-50.smat",
#     # "study-20-draft-150.smat",
# ]

# gs = [
#     # "commutes-all.smat","modmexico-city.smat", "covidflows-2020_08_31-filtered-20.smat", #row1
#     # "cn-moduillinois", "cn-Penn", "cn-modWiscon", #row2...
#     # "dblp","enron","anon", #row 3
#     # "cit-HepPh", "slashdot", "flickr", 
#     # "geometric",
#     # "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.9-0.0-100-5",
#     "study-11-2023-0-noweak.smat",
#     # "study-11-2023-1-longrange-1.smat",
#     # "study-11-2023-1-longrange-2.smat",
#     # "study-11-2023-1-longrange-3.smat",
#     # "study-11-2023-1-longrange-5.smat",
#     # "study-11-2023-1-longrange-8.smat",
#     # "study-11-2023-1-longrange-10.smat",
#     # "study-11-2023-1-longrange-12.smat",
#     # "study-11-2023-1-longrange-15.smat",
#     "study-11-2022-1.smat",
#     "study-11-2022-10.smat",
#     "study-11-2022-20.smat",
#     "study-11-2022-30.smat",
#     "study-11-2022-40.smat",
#     "study-11-2022-50.smat",
# ]

# gs = [  
#             # "modmexico-city.smat",
#             # "commutes-all.smat",
#             # "covidflows-2020_08_31.smat",
#             # "covidflows-2020_08_31-filtered-20.smat",
#             "study-20-draft-150.smat",
#             # "study-25-1.smat",
#             # "study-25-2.smat",
#             # "study-25-150.smat"
#             ]
gs = [
    "internal-shuffled-cl-louvain-0-1-study-25-1.smat",
    "internal-shuffled-cl-louvain-1-1-study-25-1.smat",
    "internal-shuffled-cl-louvain-2-1-study-25-1.smat",
    "internal-shuffled-cl-louvain-3-1-study-25-1.smat"
]


# gs = ["internal-shuffled-cl-louvain-0-$gname" for gname in gs]

# gs = getgnames("study-11-2022","input/graphs/")
# push!(gs,getgnames("noweak","input/graphs/")...)

# gs = getgnames("study-25","input/graphs/")
# gs = getgnames("study-25-150","input/graphs/")

# gs = ["study-24-2.smat","study-24-100.smat"]
# gs = map(x->getgnames(x,"input/graphs/")[1],gs)

# gs = getgnames("longrange-5-","input/graphs/")

#sort by number of edges in graph
p = sortperm(map(x->get_graph_stats(x,gpath="input/graphs/")[2],gs))
gs = gs[p]

gnames = Vector{Vector{String}}()
for g in gs
    tmp = include_graph(g,gpath)
    filter!(x->!occursin("triangle",x),tmp)
    push!(gnames,tmp)
end


ngraphs = sum(length.(gnames))

# diffusions with a fixed beta
for method in ["seir"]#,"sir"]
    for (i,g) in enumerate(gs)
        dst = joinpath(parentDst,"$(g[1:end-5])/diffusions/uniform/")
        # println("working on betas for $i of $(length(gs)) - gname: $g")
        for (j,gname) in enumerate(gnames[i])
            println("$i of $(length(gs)) - gname: $(g[1:end-5]), $j of $(length(gnames[i])) - $(gname[1:end-5])")
            curr_betas,betabool = _checkgraph(gname,bs;dst=dst,method=method,qpercents=qpercents,nnodes=nnodes)
            if length(curr_betas)>0
                qdiffusion1(nnodes,ktrials,[gname],betas=curr_betas,gpath=gpath,qpercents=qpercents,tmax=tmax,
                            dst=dst,inflb=1,nodesampletype="uniform",method=method,node_rseed=71)
            end
        end
    end
end

#=
R0s = [1.1, 5, 10, 20, 25, 26, 27, 28, 29, 30, 35, 40, 45, 50, 60, 75]


eig_data = readdlm("pipeline/data/dominant-eigvals-0.txt")
eig_data = Dict(zip(eig_data[:,1],eig_data[:,2]))

function beta_from_r0(r0,lam,gamma::Float64=5e-2)
    #r0 = lam*beta/gamma
    return r0*gamma/lam
end

#for R0 based epidemics
# for method in ["seir"]#,"sir"]
#     for (i,g) in enumerate(gs)
#         dst = joinpath(parentDst,"$(g[1:end-5])/diffusions/uniform/scratch/")
#         # println("working on betas for $i of $(length(gs)) - gname: $g")
#         for (j,gname) in enumerate(gnames[i])
#             println("$i of $(length(gs)) - gname: $g, $j of $(length(gnames[i])) - $gname")
#             #finding values of beta 
#             g_lambda = eig_data[gname]
#             bs = [beta_from_r0(x,g_lambda) for x in R0s]
#             filter!(x->x<1, bs)

#             curr_betas,betabool = _checkgraph_new(gname,bs;dst=dst,method=method,qpercents=qpercents,nnodes=nnodes)
#             if length(curr_betas[betabool])>0
#                 #for qpercent 
#                 qdiffusion1(nnodes,ktrials,[gname],betas=curr_betas[betabool],gpath=gpath,qpercents=qpercents,tmax=tmax,
#                             dst=dst,inflb=1,nodesampletype="uniform",method=method,node_rseed=71)
#                 #for explicitly defined quarantine capacities
#                 # qdiffusion1(nnodes,ktrials,[gname],betas=curr_betas[betabool],gpath=gpath,qsizes=qcaps,tmax=tmax,
#                 #             dst=dst,inflb=1,nodesampletype="uniform",method=method,node_rseed=71)
#             end
#         end
#     end
# end


#=
#specific to study graphs 
#we want to control for changes in average degree. 
#specify beta for original network (no weak ties)
#then compute betas for other networks so that beta*davg is preserved
qcaps = [0,25,50,75,100,125,150,200,250,300,400,500,1000,2000,3000,4000,5000]

#betas for base graph 
bs = unique(vcat(1e-3:1e-3:1e-2,2e-2:1e-2:1e-1))

Abase = loadGraph(gs[1],"input/graphs/")
davg_base = nnz(Abase)/lastindex(Abase,1)

for method in ["seir"]#,"sir"]
    for (i,g) in enumerate(gs)
        dst = joinpath(parentDst,"$(g[1:end-5])/diffusions/uniform/scratch/")
        # println("working on betas for $i of $(length(gs)) - gname: $g")
        #compute graph specific betas
        Atmp = loadGraph(g,"input/graphs/")
        davg_tmp = nnz(Atmp)/lastindex(Atmp,1)

        bs_tmp = bs.*(davg_base/davg_tmp)
        @show (g,bs_tmp)

        for (j,gname) in enumerate(gnames[i])
            println("$i of $(length(gs)) - gname: $g, $j of $(length(gnames[i])) - $gname")
            curr_betas,betabool = _checkgraph(gname,bs_tmp;dst=dst,method=method,qpercents=qcaps,nnodes=nnodes)
            if length(curr_betas[betabool])>0
                # qdiffusion1(nnodes,ktrials,[gname],betas=curr_betas[betabool],gpath=gpath,qpercents=qpercents,tmax=tmax,
                #             dst=dst,inflb=1,nodesampletype="uniform",method=method,node_rseed=71)
                qdiffusion1(nnodes,ktrials,[gname],betas=curr_betas[betabool],gpath=gpath,qsizes=qcaps,tmax=tmax,
                            dst=dst,inflb=1,nodesampletype="uniform",method=method,node_rseed=71)
            end
        end
    end
end
=#
=#