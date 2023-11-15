# figure for spatial graphs as we increase iterations

parent_dir = "/p/mnt/scratch/network-epi/"
include(joinpath(parent_dir,"code","make-plotting-data.jl"))

fig_dst = joinpath(parent_dir,"code","paper-figs","spatial-figs")

#generating and saving diffusion plots
#do this with julia version 1.5.1 since this throws an error for the font family on v1.8

# xpts,ypts = [12.5;13.5],[0.9;16.1]
# plot!(f,xpts,[ypts[2];ypts[2]],c=:black,leg=false)
# plot!(f,xpts,[ypts[1];ypts[1]],c=:black,leg=false)
# plot!(f,[xpts[1];xpts[1]],ypts,c=:black,leg=false)
# plot!(f,[xpts[2];xpts[2]],ypts,c=:black,leg=false)

#do base graph
gname = "study-11-2022-1.smat"
beta = 0.03

full,aggregated = qpercent_contours(gname,betas=[beta])
f = aggregated[1]
plot!(f,framestyle=:none,xlabel="",ylabel="",title="",
    margins=-20Measures.mm,colorbar=:none,
    dpi=300)

plot!(f,title="")

Plots.savefig(f,joinpath(fig_dst,"heatmap-$(gname[1:end-5]).png"))
#do other graphs 
@showprogress for steps in 5:5:50
    gname = "study-11-2022-$steps.smat"
    beta = 0.03

    full,aggregated = qpercent_contours(gname,betas=[beta])
    f = aggregated[1]
    plot!(f,framestyle=:none,xlabel="",ylabel="",title="",
        margins=-20Measures.mm,colorbar=:none,
        dpi=300)

    Plots.savefig(f,joinpath(fig_dst,"heatmap-$(gname[1:end-5]).png"))
end

#generate ncp plots

include(joinpath(parent_dir,"code","ncp","ncpplots1.jl"))
include(joinpath(parent_dir,"code","ncp","ncp-acl.jl"))
gname = "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17.smat"


function graphSize(gname)
    A = loadGraph(joinpath("pipeline/graphs/",gname))
    return size(A,1)
end


function load_ncp(gname,datadir="pipeline/data/")
    n = graphSize(gname)
    
    base_graph = canonical_graph_name(gname)
    
    ncp,headerinfo = readdlm(joinpath(datadir,base_graph[1:end-5],"ncpdata","ncpinfo-$(gname[1:end-5]).txt"),',',header=true)
    ncp = DataFrame(ncp,vec(headerinfo))
    
    x = map(x->min(ncp.size[x],n-ncp.size[x]),1:size(ncp,1))
    y = ncp.cond
  
    return x,y
end

#assemble plots in inkscape or other post-processing platform

gname = "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17.smat"
x,y = load_ncp(gname)
f = myncpplot1(x,y,ylims=(8e-4,1.2),xlims=(5,2e4),plotmin=false,plotmedian=false)
plot!(f,grid=false,dpi=200)
Plots.savefig(f,joinpath(fig_dst,"ncp-$(gname[1:end-5]).png"))


gname = "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5.smat"
x,y = load_ncp(gname)
f = myncpplot1(x,y,ylims=(8e-4,1.2),xlims=(5,2e4),plotmin=false,plotmedian=false)
plot!(f,grid=false,dpi=200)
Plots.savefig(f,joinpath(fig_dst,"ncp-$(gname[1:end-5]).png"))



A = loadGraph("study-11-2022-50.smat","input/graphs/")
nnz(A)/size(A,1)