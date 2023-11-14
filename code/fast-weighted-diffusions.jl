##
using Distributed,ProgressMeter
# addprocs(2)
#allows for execution from command line as well as an ide so files can be run modularlly
mainDir = joinpath(split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))])

@everywhere mainDir = $mainDir 
@everywhere include(joinpath(mainDir,"code/fast-diffusion1.jl"))
@everywhere include(joinpath(mainDir,"code/graph-io.jl"))
@everywhere include(joinpath(mainDir,"code/data-io.jl"))


gpath = joinpath(mainDir,"pipeline/graphs/")
parentDst = joinpath(mainDir,"pipeline/data/")


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

bs = unique(vcat(1e-4:1e-4:1e-3,1e-3:1e-3:1e-2,2e-2:1e-2:1e-1))
qpercents = collect(0:15)
nnodes,ktrials = 50,1
tmax = 20000

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


gs = ["study-20-draft-150"]
gs = vcat(map(x->getgnames(x,"input/graphs/"),gs)...)

#more params for commutes, Sparsified networks, and anon
#sort by number of edges in graph
p = sortperm(map(x->get_graph_stats(x,gpath="input/graphs/")[2],gs))
gs = gs[p]

gnames = Vector{Vector{String}}()
for g in gs
    tmp = include_graph(g,gpath)
    filter!(x->occursin("triangle",x)||(!occursin("rewired-",x)&&!occursin("er-",x)),tmp)
    push!(gnames,tmp)
end

#triangle weighted for original graph 
ngraphs = sum(length.(gnames))


## betas
for method in ["seir","sir"]
    for (i,g) in enumerate(gs)
        dst = joinpath(parentDst,"$(g[1:end-5])/diffusions/new-triangle-weighted/")
        # println("working on betas for $i of $(length(gs)) - gname: $g")
        for (j,gname) in enumerate(gnames[i])
            println("$i of $(length(gs)) - gname: $g, $j of $(length(gnames[i])) - $gname")
            curr_betas,betabool = _checkgraph(gname,bs;dst=dst,method=method,qpercents=qpercents,nnodes=nnodes)
            if length(curr_betas[betabool])>0
                triangle_qdiffusion(nnodes,ktrials,[gname],betas=curr_betas[betabool],gpath=gpath,qpercents=qpercents,tmax=tmax,
                            dst=dst,inflb=1,nodesampletype="uniform",method=method,node_rseed=71)
            end
        end
    end
end
