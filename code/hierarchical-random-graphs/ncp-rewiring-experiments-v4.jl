#= 
ncp rewiring version 4
    increase max set size to 5000 when preprocessing sets.
    basically, want to make sure we are filtering out sets 
    that have small conductance due to denominator in cut/vol.
=#


#script to perform ncp rewirings and save to scratch directory 
using Distributed 

mainDir = joinpath("/",split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))]...)
@everywhere mainDir = $mainDir

@everywhere include(joinpath(mainDir,"code/hierarchical-random-graphs/ncp-rewiring.jl"))

dstDir =  joinpath(mainDir,"scratch","ncp-rewired-graphs-v4")
@everywhere dstDir = $dstDir 

gpath = joinpath(mainDir,"input/graphs/")
@everywhere gpath = $gpath

#gnames 
gs = [
        "study-20-draft-150",
        "study-25-1",
        "study-25-2",
        "study-25-150",
        "cit-hepph",
        "email-enron",
        "filtered",
        "geometric",
        "mexico-city",
        ]

#parameters 
rewiring_ps = vec([0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.25 0.5 1.0]) #excluding 10.0 and 100.0
# rewiring_ps = vec([])

gnames = first.(getgnames.(gs,gpath))
filter!(x->endswith(x,".smat"),gnames)

params = Dict{String,Vector{Float64}}()
for gname in gnames 
    params[gname] = deepcopy(rewiring_ps)
end

#filter out graphs we've already rewired 
function filter_rewired_graphs(gname::String,gdst::String)
    graphnames = include_graph(gname,gdst)
    bool = all(map(x->"ncp-rewired-$(round(x*100,digits=2))-v4-$gname" in graphnames,params[gname]))
    return !bool 
end

filter!(x->filter_rewired_graphs(x,dstDir),gnames)
ngraphs = length(gnames)

#function for passing to pmap 
@everywhere function ncp_rewiring(
            gname::String,
            p::Float64=1.0,
            rseed::Int=-1
    )
    @show "working on $gname and $p"
    if  rseed<0
        Random.seed!(time_ns())
    else
        Random.seed!(rseed)
    end
   
    #shorter graph name 
    g = gname[1:end-5]

    #load graph 
    A = loadGraph(gname,"input/graphs/")

    #load precomputed NCP and sets 
    ncploc = joinpath(mainDir,"pipeline/data/$g/ncpdata/")
    ncp,headerinfo = readdlm(ncploc*"ncpinfo-$g.txt",',',header=true)
    ncp = DataFrame(ncp,vec(headerinfo))

    ncpmat = MatrixNetworks.readSMAT(ncploc*"ncpmat-$g.smat")
    sets = ncpmat_to_sets(ncpmat)
    
    #preprocess sets.  
    ncp,sets = preprocess_sets!(ncp,sets,A,5000)

    #ensure that sets touch all the nodes 
    s = Set{Int}()
    for x in sets
        union!(s,x)
    end
    @assert(length(s)==lastindex(A,1),"Sets don't cover every node, some edges are static and won't be modified.")
    
    #perform rewiring 
    nsteps = round(Int,p*nnz(A))
    new_edges = conductance_rewiring(A,ncp,sets,nsteps)
    X = edges2graph(new_edges)

    #write to dst 
    fname = joinpath(dstDir,"ncp-rewired-$(round(p*100,digits=2))-v4-$gname")  
    writeSMAT(X,fname)

    #warnings
    if nnz(X)!=nnz(A)
        @warn("different number of edges in original and rewired graphs")
        @show nnz(A)
        @show nnz(X)
        @show fname
    end
end

#main call 
parameters = vcat(Iterators.product(gnames,rewiring_ps)...)
@showprogress pmap(p->ncp_rewiring(p[1],p[2]),parameters)