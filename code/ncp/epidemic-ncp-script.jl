using Distributed, ProgressMeter
# addprocs(50)

parent_dir = "/p/mnt/scratch/network-epi/"
println("Loading scripts on workers")
@everywhere include(joinpath($parent_dir,"code/fast-diffusion.jl"))
@everywhere include(joinpath($parent_dir,"code/graph-io.jl")) 
@everywhere include(joinpath($parent_dir,"code/data-io.jl")) 

## plotting fcns
include(joinpath(parent_dir,"code/ncp/ncpplots1.jl"))
# addprocs(15)
ENV["GKSwstype"] = "100"
gpath = "pipeline/graphs/"

@everywhere include(joinpath($parent_dir,"code/ncp/parallel-epidemic-ncp.jl"))
println("scripts loaded")

gs = readdir(joinpath(parent_dir,"input/graphs/"))
filter!(x->endswith(x,".smat"),gs)
#ignoring these for now 
filter!(c->!occursin("livejournal",lowercase(c)),gs)

#sort by number of edges in graph
p = sortperm(map(x->get_graph_stats(x,gpath="input/graphs/")[2],gs))

#TODO filter out networks whose epidemic ncp we already have 


# gs  = ["cit-HepPh.smat","cit-Patents.smat"]

# gs = ["rewired-10000.0-modmexico-city.smat",
#     "rewired-10000.0-cn-modUIllinois20.smat",
#     "er-10000.0-dblp-cc.smat"]

total_trials = 50000
# gs = ["study-11-2023-1-longrange-1.smat",
# "study-11-2023-1-longrange-2.smat",
# "study-11-2023-1-longrange-3.smat",
# "study-11-2023-1-longrange-5.smat",
# "study-11-2023-1-longrange-8.smat",
# "study-11-2023-1-longrange-10.smat",
# "study-11-2023-1-longrange-12.smat",
# "study-11-2023-1-longrange-15.smat",
# ]
gs = getgnames("study-24","input/graphs/")
push!(gs,getgnames("study-25","input/graphs/")...)
# gs = [
#     "study-25-1.smat",
#     "study-25-2.smat",
#     "study-25-5.smat",
#     "study-25-10.smat",
#     "study-25-20.smat",
#     "study-25-25.smat",
#     "study-25-50.smat",
#     "study-25-100.smat",
#     "study-25-150.smat",
#     "study-24-1.smat",
#     "study-24-2.smat",
#     "study-24-5.smat",
#     "study-24-10.smat",
#     "study-24-20.smat",
#     "study-24-30.smat",
#     "study-24-50.smat",
#     "study-24-70.smat",
#     "study-24-100.smat"
# ]

for model in ["seir"]#,"sir"]
    for gname in gs
        g = canonical_graph_name(gname)
        dst = joinpath(parent_dir,"pipeline/data/$(g[1:end-5])/ncpdata/")
        fnames = readdir(dst)
        #check number of nodes to decide if we want to do extremal graphs too
        basegraph_path = joinpath(parent_dir,"input","graphs",gname)
        # if parse(Int,split(readline(basegraph_path)," ")[1])<1000000 #more than 1M nodes
        #     gnames = [gname, "rewired-10000.0-$gname", "er-10000.0-$gname"]
        # else 
        #     gnames = [gname]
        # end
        # gnames = [gname]
        gnames = [gname, "rewired-10000.0-$gname", "er-10000.0-$gname"]

        for h in gnames
            println("working on $(uppercase(model)) epidemic ncp for $h")
            if "ncpinfo-epidemic-subsampling-4-$model-$total_trials-$(h[1:end-5]).txt" in fnames
                println("epidemic ncp data already exists for $h")
            else
                new_parallel_ncp_code3(h,total_trials,model=model,dst=dst)
            end
            # new_parallel_ncp_code3(h,total_trials,model=model,dst=dst)
        end
    end
end

#saving images 
# total_trials = 50000
# for model in ["seir","sir"]
#     for gname in gs
#         dst = joinpath(parent_dir,"pipeline/data/$(gname[1:end-5])/ncpdata/")
#         imgdst = joinpath(parent_dir,"pipeline/data/$(gname[1:end-5])/imgs/")
#         fnames = readdir(dst)
#         for h in [gname, "rewired-10000.0-$gname", "er-10000.0-$gname"]
#             println("working on $(uppercase(model)) for $h")
#             fname = "ncpinfo-epidemic-subsampling-$model-$total_trials-$(h[1:end-5]).txt"
#             #load data
#             A = loadGraph(h,"pipeline/graphs/")
#             ncp,headerinfo = readdlm(joinpath(dst,fname),',',header=true)
#             ncp = DataFrame(ncp,vec(headerinfo))
#             x,y = map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1)),ncp.cond
#             f = myncpplot1(x,y,
#                 plotmin=false,plotmedian=true,
#                 ylims=(1e-4,1.1),
#                 xlims=(5,maximum(x)*5))
                
#             title!(f,"$(h[1:end-5]) Epidemic NCP - Subsampling")
#             xlabel!(f,"Size")
#             ylabel!(f,"Conductance")
#             savefig(f,joinpath(imgdst,"hexbinncp-epidemic-subsampling-$model-$total_trials-$(h[1:end-5]).png")) 
#         end
#     end
# end


#=
#custome calculation for small spatial network  
using Measures
gname = getgnames("lfr","input/graphs/")[end]
A = loadGraph("input/graphs/$gname")
dst = joinpath("pipeline/data/$(gname[1:end-5])/ncpdata/")

ncp,headerinfo = readdlm(joinpath(dst,"ncpinfo-epidemic-subsampling-4-seir-50000-$(gname[1:end-5]).txt"),',',header=true)
# ncp,headerinfo = readdlm(joinpath(dst,"ncpinfo-epidemic-subsampling-4-seir-2000-$(gname[1:end-5]).txt"),',',header=true)
ncp = DataFrame(ncp,vec(headerinfo))

# f = myncpplot1(map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1)),ncp.cond,
#     plotmin=false,plotmedian=false,nbins=(80,100),
#     xlims=(5,800),
#     ylims=(1e-2,1.1))

# xs = xlims(f)
# ys = ylims(f)
# f = plot([],[],xlims=xs,ylims=ys,leg=false,
#     xscale=:log10,yscale=:log10)
using PerceptualColourMaps
f = plot()
myhexbin_conditional!(f,map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1)),ncp.cond;
        ylims=(1e-4,1.1),nbins=(80,80),
        color=cgrad(cmap("D06")),colorbar=true)

f = plot()
myhexbin_conditional!(f,map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1)),ncp.cond;
    ylims=(1e-4,1.1),nbins=(80,80),
    color=cgrad([:black,:blue, :orange],[0.1,0.3,0.5,0.7,0.9]),
    colorbar=true)


# f = myhexbin_conditional(map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1)),ncp.cond,
#     nbins=(80,80),
#     ylims=(1e-2,1.1))

# f

f = Plots.plot(10.0.^(res[1]),10.0.^(res[2]),fill_z=log10.(res[4]),linecolor=nothing,
        seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=false,
        xlims=(2.8,100000),ylims=(1e-3,1.5))

gnames = ["dblp","enron",
        "anon",
        "mexico","commutes","covid",
        "study-20-draft-150.smat","study-11-2022-45",
        "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5",
        "cn-modUill"]

gnames = vcat(getgnames.(gnames,"input/graphs/")...)


figs = []
for gname in gnames
    println("working on $gname")
    A = loadGraph("input/graphs/$gname")
    dst = joinpath("pipeline/data/$(gname[1:end-5])/ncpdata/")
    ncp,headerinfo = readdlm(joinpath(dst,"ncpinfo-epidemic-subsampling-seir-50000-$(gname[1:end-5]).txt"),',',header=true)
    ncp = DataFrame(ncp,vec(headerinfo))

    f = myhexbin_conditional(map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1)),ncp.cond,
    nbins=(80,100),ylims=(1e-3,1.2))
    Plots.plot!(f,title="$(gname[1:end-5])")

    push!(figs,f)
end

figs[1]
figs[2]
figs[3]
figs[4]
figs[5]
figs[6]
figs[7]
figs[8]
figs[9]
figs[10]
figs[11]


# using Measures
# gname = getgnames("slash","input/graphs/")[1]
# A = loadGraph(gname,"input/graphs/")
# fdst = "/p/mnt/scratch/network-epi/scratch/figures/"
# #epidemic ncp 
# dst = joinpath("pipeline/data/$(gname[1:end-5])/ncpdata/")
# ncp,headerinfo = readdlm(joinpath(dst,"ncpinfo-epidemic-seir-60000-$(gname[1:end-5]).txt"),',',header=true)
# ncp = DataFrame(ncp,vec(headerinfo))
# f = myncpplot1(map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1)),ncp.cond,
#     plotmin=false,plotmedian=true,ylims=(1e-3,1.1))

# title!(f,"$(gname[1:end-5]) Epidemic NCP")
# xlabel!(f,"Size")
# ylabel!(f,"Conductance")
# xlims!(f,(6, 50000))

# #subsampling
# ncp,headerinfo = readdlm(joinpath(dst,"ncpinfo-epidemic-subsampling-seir-50000-$(gname[1:end-5]).txt"),',',header=true)
# ncp = DataFrame(ncp,vec(headerinfo))
# f = myncpplot1(map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1)),ncp.cond,
#     plotmin=false,plotmedian=true,ylims=(1e-3,1.1))
    
# title!(f,"$(gname[1:end-5]) Epidemic NCP - Subsampling")
# xlabel!(f,"Size")
# ylabel!(f,"Conductance")
# xlims!(f,(6, 50000))
# # savefig(f,joinpath(fdst,"hexbinncp-epidemic-seir-60000-$(gname[1:end-5]).png"))


# #acl ncp 
# ncp,headerinfo = readdlm(joinpath(dst,"ncpinfo-$(gname[1:end-5]).txt"),',',header=true)
# ncp = DataFrame(ncp,vec(headerinfo))
# inds = findall(ncp.cond.<0)

# ncp = ncp[setdiff(1:size(ncp,1),inds),:]
# f = myncpplot(map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1)),ncp.cond,
#     plotmin=false,plotmedian=true)

# title!(f,"Flickr ACL NCP")
# xlabel!(f,"Size")
# ylabel!(f,"Conductance")
# ylims!(f,(0.0007600059937660782, 1.3229577599736124))
# xlims!(f,(6, 500000))

# savefig(f,joinpath(fdst,"hexbinncp-acl-$(gname[1:end-5]).png"))

# f

# pwd()
# include("../../code/graph-io.jl")


# getgnames("covid","input/graphs/")

# gnames = ["loc-gowalla.smat"; "loc-brightkite.smat"; "covidflows-2020_08_31.smat";
#             "covidflows-2020_08_31-filtered-20.smat"; "soc-Epinions1.smat"]

# gname = gnames[1]

# #ACL NCP
# parentDir = "/p/mnt/scratch/network-epi/pipeline/data/"
# joinpath(parentDir,"$(gname[1:end-5])","ncpdata")


dst = joinpath(parent_dir,"pipeline/data/$(gs[1][1:end-5])/ncpdata/")
fnames = readdir(dst)
filter!(c->occursin("worker",c),fnames)
=#