using Plots
using Random
using ProgressMeter

mainDir = joinpath("/",split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))]...)
dstDir =  joinpath(mainDir,"input","graphs")

include(joinpath(mainDir,"code","graph-io.jl"))
include(joinpath(mainDir,"code","hierarchical-random-graphs/nested-sbm.jl"))
# include("../ncp/ncp-acl.jl")

#generate and save a few networks w/ ~50k nodes 
groupsize = [50]
nlayers = [10] 
avgd = [20] #should be LESS THAN groupsize or we truncate to min(groupsize-1,avgd)
decay_factor=[0.7, 0.5, 0.4, 0.2]

#set random seeds for networks 
parameters = Iterators.product(groupsize,nlayers,avgd,decay_factor)

rseeds = 1:length(parameters)
@showprogress for (i,(gs,nl,d,decay)) in enumerate(parameters)
    @show (gs,nl,d,decay)
    Random.seed!(rseeds[i])
    A = nested_sbm(gs,nl,d,decay)
    A = largest_component(A)[1]

    #save graph 
    fname = "nested-sbm-$gs-$nl-$d-$decay-$(rseeds[i]).smat"
    writeSMAT(A,joinpath(dstDir,fname))    
end


#=
ncp,sets = make_ncp(A)
x = map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1))
y = ncp.cond

extrema(x)
extrema(y)
# #from the "missed-sets-v2.jl" script
# nlines = 2
# f = make_base_plot(nlines,(0.8,1e4),10)
# plot!(f,size=(400,400))
# ybounds = (10.0^(-nlines-1),2.0)
# myhexbin2!(f,x,y)
# plot!(f;ylims=ybounds,xlims=(0.8,1e4),clims=(0,1))

f = scatter(x,y,xscale=:log10,yscale=:log10,leg=false,
    xlabel="Size",
    ylabel="Conductance",
    guidefontsize=16,
    tickfontsize=12,
    xlims=(0.8,maximum(x)*2),
    ylims=(minimum(y)/8,1.05),
    top_margin=4Measures.mm)

title_str = "Nested SBM PPR-NCP"
title_str*="\nnnodes:$(lastindex(A,1)), avgd=$(round(nnz(A)/lastindex(A,1),digits=2))"
title_str*="\ngroupsize: $groupsize, nlayers: $nlayers, decay rate: $decay_factor"
title!(f,title_str)
=#