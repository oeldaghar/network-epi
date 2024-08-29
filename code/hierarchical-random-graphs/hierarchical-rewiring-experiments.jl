using Plots
using Measures
using Random
using ProgressMeter

mainDir = joinpath("/",split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))]...)
dstDir =  joinpath(mainDir,"input","graphs")

include(joinpath(mainDir,"code","graph-io.jl"))
include(joinpath(mainDir,"code","ncp","ncp-acl.jl"))
include(joinpath(mainDir,"code","hierarchical-rewiring","ncp-rewiring.jl"))
# include(joinpath(mainDir,"code","hierarchical-random-graphs/hierarchical-rewiring.jl"))


#parameters for graph experiments 
