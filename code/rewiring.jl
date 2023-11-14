# script for rewiring graphs. see https://arxiv.org/abs/1202.3473 for
#discussion about threshold for rewiring.

#allows for execution from command line as well as an ide so files can be run modularlly
mainDir = joinpath(split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))])

include(joinpath(mainDir,"code","graph-io.jl"))
include(joinpath(mainDir,"code","rewiring-functions.jl"))

using Distributions, Random
gpath = joinpath(mainDir,"input/graphs/")
gdst = joinpath(mainDir,"pipeline/graphs/")

rewiring_ps = vec([0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.25 0.5 1.0 10.0 100.0])
#for large graphs (flickr)
alt_rewiring_ps = vec([0.001 0.01 0.1 1.0 100.0]) 


gnames = readdir(gpath)
filter!(x->endswith(x,".smat"),gnames)
params = Dict{String,Vector{Float64}}()
for gname in gnames 
    params[gname] = deepcopy(rewiring_ps)
end

#special handling for large graphs (flickr)
params[getgnames("flickr",gpath)[1]] = deepcopy(alt_rewiring_ps)

#filter out graphs we've already rewired 
function filter_rewired_graphs(gname::String,gdst::String)
    graphnames = include_graph(gname,gdst)
    bool = all(map(x->"rewired-$(round(x*100,digits=2))-$gname" in graphnames,params[gname]))
    bool = bool && all(map(x->"er-$(round(x*100,digits=2))-$gname" in graphnames,params[gname]))
    return !bool 
end

filter!(x->filter_rewired_graphs(x,gdst),gnames)
ngraphs = length(gnames)

#degree preserving and ER rewirings 
c = [1]
for (i,gname) in enumerate(gnames)
    rewiring_params = params[gname]
    #proportion of edges rewired. 1.0 = |E| sampled edges  
    nrewirings = 2*length(rewiring_params)
    Random.seed!(7)
    seeds = abs.(rand(Int,nrewirings))

    println("working on $i of $ngraphs : $gname")
    A = loadGraph(gname,gpath)
    rewire_steps = nnz(A)#2*ceil(Int,nnz(A)*log(nnz(A)))

    #reset seeds
    c[1] = 1

    #rewired graphs
    println("working on degree preserving rewirings")
    @showprogress for p in rewiring_params
        Random.seed!(seeds[c[1]])
        c[1] += 1
        B = rewire_graph(A,ceil(Int,p*rewire_steps))
        #save graph
        writeSMAT(B,gdst*"rewired-$(round(p*100,digits=2))-$gname")
    end

    println("working on ER rewiring")
    @showprogress for p in rewiring_params
        Random.seed!(seeds[c[1]])
        c[1] += 1
        B = er_rewire_graph(A,ceil(Int,p*rewire_steps))
        #save graph
        writeSMAT(B,gdst*"er-$(round(p*100,digits=2))-$gname")
    end
end

#triangle rewirings 
total_steps = 1000000
rseed = 7

gnames = readdir(gpath)
gnames = getgnames("study-25-",gpath)[end:end]
filter!(x->endswith(x,".smat"),gnames)
for gname in gnames 
    params[gname]=collect(0.1:0.1:1.0) 
end

#custom params for certain graphs 
params[getgnames("flickr",gpath)[1]] = [0.1;0.5;1.0]

#filter out graphs we've already rewired
function filter_rewired_graphs(gname::String,gdst::String)
    graphnames = include_graph(gname,gdst)
    bool = all(map(x->"triangle-rewired-$x-$rseed-$gname" in graphnames,params[gname]))
    return !bool 
end

filter!(x->filter_rewired_graphs(x,gdst),gnames)
ngraphs = length(gnames)

for (i,gname) in enumerate(gnames) 
    ps = params[gname]
    Random.seed!(rseed)
    seeds = abs.(rand(Int,length(ps)))
    
    println("working on triangle rewiring for graph $i of $ngraphs : $gname")
    A = loadGraph(gname,gpath)

    @showprogress for (i,p) in enumerate(ps)
        Random.seed!(seeds[i])
        B = swap_triangles(A,Int(p*total_steps))
        writeSMAT(B,gdst*"triangle-rewired-$p-$rseed-"*gname)
    end
end
