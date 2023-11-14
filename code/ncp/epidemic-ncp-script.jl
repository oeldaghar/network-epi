using Distributed, ProgressMeter
# addprocs(50)

#allows for execution from command line as well as an ide so files can be run modularlly
mainDir = joinpath(split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))])

println("Loading scripts on workers")
@everywhere include(joinpath($mainDir,"code/fast-diffusion.jl"))
@everywhere include(joinpath($mainDir,"code/graph-io.jl")) 
@everywhere include(joinpath($mainDir,"code/data-io.jl")) 

## plotting fcns
include(joinpath(mainDir,"code/ncp/ncpplots1.jl"))
# addprocs(15)
ENV["GKSwstype"] = "100"
gpath = "pipeline/graphs/"

@everywhere include(joinpath($mainDir,"code/ncp/parallel-epidemic-ncp.jl"))
println("scripts loaded")

gs = readdir(joinpath(mainDir,"input/graphs/"))
filter!(x->endswith(x,".smat"),gs)
#ignoring these for now 
filter!(c->!occursin("livejournal",lowercase(c)),gs)

#sort by number of edges in graph
p = sortperm(map(x->get_graph_stats(x,gpath="input/graphs/")[2],gs))

#TODO filter out networks whose epidemic ncp we already have 

total_trials = 50000
gs = getgnames("study-20-draft-150","input/graphs/")

for model in ["seir"]#,"sir"]
    for gname in gs
        g = canonical_graph_name(gname)
        dst = joinpath(mainDir,"pipeline/data/$(g[1:end-5])/ncpdata/")
        fnames = readdir(dst)
    
        basegraph_path = joinpath(mainDir,"input","graphs",gname)
    
        gnames = [gname, "rewired-10000.0-$gname", "er-10000.0-$gname"]

        for h in gnames
            println("working on $(uppercase(model)) epidemic ncp for $h")
            if "ncpinfo-epidemic-subsampling-4-$model-$total_trials-$(h[1:end-5]).txt" in fnames
                println("epidemic ncp data already exists for $h")
            else
                new_parallel_ncp_code3(h,total_trials,model=model,dst=dst)
            end
        end
    end
end
