#a (hopefully) organized place to maintain plots and various fcns for plotting

using MatrixNetworks, SparseArrays, LinearAlgebra, Distributions, DelimitedFiles, ProgressMeter, Random, Plots

function read_inf_data(gname;dloc=dloc,beta=0.1,gamma=0.05,method="sir",dtype="cinfs",exp=5.0)
    if endswith(gname,".smat")
        g = gname[1:end-5]
    else
        g = gname
    end
    if beta!=0
        fname = dloc*"$dtype-$g-$method-$beta-$gamma.txt"
    else
        fname = dloc*"$dtype-$g-$method-$exp.txt"
    end
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

using StatsBase
#TODO this fcn appears to have some trouble when there is a concentration of measure for xs
# try visualizing a histogram of total infections for mexico city for a single beta to see example
# basically need a way to increase the size of the bins based on input xs
function logscale_histogram(xs,ymin=1/length(xs),normalize_y::Bool=true;logx::Bool=false,maxbins::Int=5000)
    maxbins = min(maxbins,ceil(1+max(abs.(extrema(xs))...)))
    h = float(StatsBase.fit(Histogram,xs,nbins=maxbins))
    if normalize_y
        h.weights./=length(xs)
        ymax = 1
    else
        ymax = length(xs)
    end

    #a little padding to make plots look better
    bins = [collect(h.edges[1]);maximum(h.edges[1])+step(h.edges[1])]
    w = [h.weights;0] .+ 0.5*ymin
    if logx
        f = plot(bins,w,seriestype=:barbins,label="",
                xlims=(1,2*maximum(xs)),xscale=:log10,
                ylims=(ymin,ymax),yscale=:log10,bar_width=step(h.edges[1]),lims=:round)
    else
        f = plot(bins,w,seriestype=:barbins,label="",
                xlims=(0,maximum(xs)+2),ylims=(ymin,ymax),yscale=:log10,lims=:round)
    end
    return f
end

## basic graph statistics
# deg dist, kcore dist, ncps
function degdist(gname::String,gpath::String)
    A = loadGraph(gname,gpath)
    d = vec(sum(A;dims=2))
    f = logscale_histogram(d,logx=true)
    plot!(f,xlabel="degree",ylabel="normalized frequency",
        title="$(gname[1:end-5]) - Degree Distribution")
    return f
end

function kcoredist(gname::String,gpath::String)
    A = loadGraph(gname,gpath)
    kcore = corenums(A)[1]
    f = logscale_histogram(kcore)
    plot!(f,xlabel="core number",ylabel="normalized frequency",
        title="$(gname[1:end-5]) - Kcore Distribution")
    return f
end

# gpath = "graphs/tmp/"
# gs = ["astro","enron","penn","ucf","unc"]
# for h in gs
#     dfigs = Vector()
#     kfigs = Vector()
#     gnames = getgnames(h,gpath)
#     for g in gnames
#         f = degdist(g,gpath)
#         push!(dfigs,f)
#         f = kcoredist(g,gpath)
#         push!(kfigs,f)
#     end
#     #normalize plots to have same limits
#     xmin,xmax = minimum(xlims.(dfigs)[1]),maximum(xlims.(dfigs)[2])
#     ymin,ymax = minimum(ylims.(dfigs)[1]),maximum(ylims.(dfigs)[2])
#     for (i,f) in enumerate(dfigs)
#         xlims!(f,(xmin,xmax))
#         ylims!(f,(ymin,ymax))
#         Plots.savefig(f,"purgatory/degdist-$(gnames[i][1:end-5]).png")
#     end
#
#     xmin,xmax = minimum(xlims.(kfigs)[1]),maximum(xlims.(kfigs)[2])
#     ymin,ymax = minimum(ylims.(kfigs)[1]),maximum(ylims.(kfigs)[2])
#     for (i,f) in enumerate(kfigs)
#         xlims!(f,(xmin,xmax))
#         ylims!(f,(ymin,ymax))
#         Plots.savefig(f,"purgatory/kcoredist-$(gnames[i][1:end-5]).png")
#     end
# end

## rewiring plots and heatmaps
# include("fast-diffusion.jl")
# using Plots
# using Measures
# function graph_heatmap(gname::String,beta::Float64,gamma::Float64=0.05;
#                         gpath::String="",dtype::String="tinfs",dloc::String="",
#                         method::String="sir",
#                         rps::Vector{Float64}=vec([0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.25 0.5 1.0 10.0]))
#
#     gr()
#     # do stuff
#     A = loadGraph(gname,gpath)
#     nnodes = size(A,1)
#     ths = nnodes.*collect(0.01:0.01:1.0)
#     g = gname[1:end-5]
#
#     data = Vector{Vector{Vector{Int}}}()
#     gnames = Vector{String}()
#     for p in 100 .*reverse(rps)# 0.1 0.25 0.5 1.0]# 10.0]#[1 5 10 25 50 100 1000]
#         push!(data,cumsum.(read_inf_data("rewired-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype)))
#         push!(gnames,"rewired-$p-$g")
#     end
#     push!(data,cumsum.(read_inf_data(g,dloc=dloc,beta=beta,gamma=gamma,dtype=dtype)))
#     push!(gnames,g)
#     for p in 100 .*rps# 0.1 0.25 0.5 1.0]# 10.0]#[1 5 10 25 50 100 1000]
#         push!(data,cumsum.(read_inf_data("er-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype)))
#         push!(gnames,"er-$p-rewired-$g")
#     end
#
#     #build heatmap
#     #will do color as pr( >=k% of nodes)
#     pdata = zeros(length(ths),length(gnames))
#     for graphind=1:length(gnames)
#         if dtype == "tinfs"
#             ps = map(x->findlast(x.>ths),last.(data[graphind][1:500]))
#         elseif dtype == "cinfs"
#             ps = map(x->findlast(x.>ths),maximum.(data[graphind][1:500]))
#         end
#         ps = Vector{Union{Int,Nothing}}(ps)
#         ps[ps.==nothing].= 0
#         ps = Int.(ps)
#         pdata[:,graphind] .= fit(Histogram,ps,1:length(ths)+1).weights
#         pdata[:,graphind] .= cumsum(pdata[end:-1:1,graphind])[end:-1:1]./500 #500 was the sample size for the data
#     end
#     if dtype == "tinfs"
#         ylab = "Pr(Total Infs ≥ k% of Nodes)"
#     elseif dtype == "cinfs"
#         ylab = "Pr(Maximum Infs ≥ k% of Nodes)"
#     end
#     f1 = heatmap(pdata,ylabel="Percent of Total Nodes",colorbar_title=ylab,
#                 title="$(gname[1:end-5]) - SIR($beta,$gamma)",clims=(0,1),
#                 xlabel="     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent",
#                 xrotation = 90)
#
#     if dtype == "tinfs"
#         clabel = "Fraction of Infected Nodes"
#         data = hcat(map(x->last.(x),data)...)./nnodes
#     elseif dtype == "cinfs"
#         clabel = "Normalized Maximum Infection Size"
#         data = hcat(map(x->maximum.(x),data)...)./nnodes
#     end
#     f = heatmap(data,ylabel="quarantine percentage",
#         xlabel="     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent",
#         title="$(gname[1:end-5]) - SIR($beta,$gamma)\n500 node sample",
#         colorbar_title=clabel,c=:heat,clims=(0,1),
#         yticks=(collect(250:500:8000),0:15),xrotation = 90)
#
#     plot!(f,xticks=(1:2*length(rps)+1, string.(round.(100 .*vec([reverse(rps)... 0 rps...]),digits=2))),margins=9Measures.mm)
#     plot!(f1,xticks=(1:2*length(rps)+1, string.(round.(100 .*vec([reverse(rps)... 0 rps...]),digits=2))),margins=9Measures.mm)
#
#     # plot!(f,xticks=(1:2*length(rps), string.(round.(100 .*vec([reverse(rps)... rps...]),digits=2))),margins=9mm)
#     # plot!(f1,xticks=(1:2*length(rps), string.(round.(100 .*vec([reverse(rps)... rps...]),digits=2))),margins=9mm)
#     return f,f1
# end

# function graph_beta_heatmap()
# end
# function graph_qpercent_heatmap()
# end
# ps = vec([0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.25 0.5 1.0 10.0 100.0])
#
# gpath = "graphs/real-graphs/"
# dloc = "diffusion-data/uf21/"
#
# bs = sort([collect(0.001:0.001:0.01);collect(0.02:0.01:0.1)])
#
# gs = ["cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.9-0.0-100-5",
#         "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5",
#         "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17"]
# gpath = "graphs/tmpgraphs3/"
# for h in gs
#     gnames = getgnames(h,gpath)
#     filter!(x->!startswith(x,"er-") && !startswith(x,"rewired-") && !startswith(x,"web-") && !startswith(x,"kcore"),gnames)
#     # filter!(x->!startswith(x,"cn-") && !startswith(x,"cl-") ,gnames)
#     for gname in gnames
#         g = gname[1:end-5]
#         # println("working on $gname")
#         println("working on $gname")
#         for beta in bs
#             println("working on beta: $beta")
#             #total infs plots
#             f,f1 = graph_heatmap(gname,beta,0.05,gpath=gpath,dloc="diffusion-data/$g/",dtype="tinfs",rps = ps)
#             Plots.savefig(f,"purgatory/totalinfs-heatmaps/totalinfs-$g-$beta-1.png")
#             Plots.savefig(f1,"purgatory/totalinfs-heatmaps/totalinfs-$g-$beta-2.png")
#             #maximum infs plots
#             # f,f1 = graph_heatmap(gname,beta,0.05,gpath=gpath,dloc=dloc,dtype="cinfs",rps = ps)
#             # Plots.savefig(f,"purgatory/peakinfs-$(gname[1:end-5])-$beta-1.png")
#             # Plots.savefig(f1,"purgatory/peakinfs-$(gname[1:end-5])-$beta-2.png")
#         end
#     end
# end


# A = loadGraph(getgnames("journal","graphs/all-graphs/")[1],"graphs/all-graphs/")
# gnames = getgnames("cn-penn","graphs/all-graphs/")

# beta|
# beta|         heatmap
# beta|         heatmap
# beta|
####### #rewiring #rewiring #rewiring

# gpath = "graphs/all-graphs/"
# dloc = "diffusion-data/tmp/"
# gnames = getgnames(gs[1],gpath)
# gname = gnames[1]
#
# #TODO modify this to that we get an image for each quarantine percentage
# function graph_heatmap0(gname::String,beta::Float64,gamma::Float64=0.05;
#                         gpath::String="",dtype::String="tinfs",dloc::String="",
#                         method::String="sir",
#                         rps::Vector{Float64}=vec([0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.25 0.5 1.0 10.0]))
#
#     gr()
#     A = loadGraph(gname,gpath)
#     nnodes = size(A,1)
#     ths = nnodes.*collect(0.01:0.01:1.0)
#     g = gname[1:end-5]
#
#     #horizontal - rewiring graph
#     #vertical - betas
#     data = Array{Float64,2}(undef,0,25)
#
#     gnames = Vector{String}()
#     g = gname[1:end-5]
#     # beta = 0.001
#     for beta in sort(bs)
#         tmpdata = Vector{Vector{Vector{Int}}}()
#         for p in 100 .*reverse(rps)# 0.1 0.25 0.5 1.0]# 10.0]#[1 5 10 25 50 100 1000]
#             push!(tmpdata,cumsum.(read_inf_data("rewired-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype)[1:500]))
#             push!(gnames,"rewired-$p-$g")
#         end
#         push!(tmpdata,cumsum.(read_inf_data(g,dloc=dloc,beta=beta,gamma=gamma,dtype=dtype)[1:500]))
#         push!(gnames,g)
#         for p in 100 .*ps# 0.1 0.25 0.5 1.0]# 10.0]#[1 5 10 25 50 100 1000]
#             push!(tmpdata,cumsum.(read_inf_data("er-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype)[1:500]))
#             push!(gnames,"er-$p-rewired-$g")
#         end
#
#         data = vcat(data,hcat(map(x->last.(x),tmpdata)...)./nnodes)
#     end
#
#     clabel = "Fraction of Infected Nodes"
#
#     f = heatmap(data,ylabel="",
#         xlabel="     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent",
#         title="$(gname[1:end-5]) - SIR(β,$gamma)\n500 node sample",
#         colorbar_title=clabel,c=:heat,clims=(0,1),
#         yticks=(250:500:19*500+250,sort(bs)),xrotation=90)
#
#     plot!(f,xticks=(1:2*length(ps)+1, string.(round.(100 .*vec([reverse(ps)... 0 ps...]),digits=2))),margins=9Measures.mm)
#
#     return f
# end
#
# #zoom in on case with no quarantining
# function graph_heatmap_noq(gname::String,betas::Vector{Float64},gamma::Float64=0.05;
#                         gpath::String="",dtype::String="tinfs",dloc::String="",
#                         method::String="sir",
#                         rps::Vector{Float64}=vec([0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.25 0.5 1.0 10.0]))
#
#     gr()
#     A = loadGraph(gname,gpath)
#     nnodes = size(A,1)
#     ths = nnodes.*collect(0.01:0.01:1.0)
#     g = gname[1:end-5]
#
#     #horizontal - rewiring graph
#     #vertical - betas
#     data = Vector{Vector{Vector{Int}}}()
#     tmpdata = Vector{Vector{Int}}()
#     gnames = Vector{String}()
#     for p in 100 .*reverse(rps)# 0.1 0.25 0.5 1.0]# 10.0]#[1 5 10 25 50 100 1000]
#         push!(data,cumsum.(read_inf_data("rewired-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype)[1:500]))
#         push!(gnames,"rewired-$p-$g")
#     end
#     push!(data,cumsum.(read_inf_data(g,dloc=dloc,beta=beta,gamma=gamma,dtype=dtype)[1:500]))
#     push!(gnames,g)
#     for p in 100 .*rps# 0.1 0.25 0.5 1.0]# 10.0]#[1 5 10 25 50 100 1000]
#         push!(data,cumsum.(read_inf_data("er-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype)[1:500]))
#         push!(gnames,"er-$p-rewired-$g")
#     end
#
#     if dtype == "tinfs"
#         clabel = "Fraction of Infected Nodes"
#         data = hcat(map(x->last.(x),data)...)./nnodes
#     elseif dtype == "cinfs"
#         clabel = "Normalized Maximum Infection Size"
#         data = hcat(map(x->maximum.(x),data)...)./nnodes
#     end
#     f = heatmap(data,ylabel="quarantine percentage",
#         xlabel="     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent",
#         title="$(gname[1:end-5]) - SIR($beta,$gamma)\n500 node sample",
#         colorbar_title=clabel,c=:heat,clims=(0,1),
#         yticks=(collect(250:500:8000),0:15),xrotation = 90)
#         # xticks=(collect(1:2*length(ps)+1), string.(vec([reverse(ps)... 0 ps...]))))#clims=(0,1))
#     # plot!(f,yticks=(collect(250:500:8000),0:15))
#     # plot!(f,xticks=([1;5;9],["Original";"Rewired";"ER"]))
#     # plot!(f,xticks=(1:2*length(ps)+1, string.(vec([reverse(ps)... 0 ps...]))))
#
#     plot!(f,xticks=(1:2*length(rps)+1, string.(round.(100 .*vec([reverse(rps)... 0 rps...]),digits=2))),margins=9Measures.mm)
#     plot!(f1,xticks=(1:2*length(rps)+1, string.(round.(100 .*vec([reverse(rps)... 0 rps...]),digits=2))),margins=9Measures.mm)
#
#     # plot!(f,xticks=(1:2*length(rps), string.(round.(100 .*vec([reverse(rps)... rps...]),digits=2))),margins=9mm)
#     # plot!(f1,xticks=(1:2*length(rps), string.(round.(100 .*vec([reverse(rps)... rps...]),digits=2))),margins=9mm)
#     return f,f1
# end
#
# function peakInfsCompPlot(gname::String,timescale::Int;gpath::String="",dloc::String="",beta::Float64=0.05,gamma::Float64=0.05,method::String="sir",dst::String="",nnodes::Int=1)
#     gnames = ["rewired-1000.0-$gname";gname;"er-1000.0-$gname"]
#     ys = Vector{Vector{Float64}}()
#     for g in gnames
#         push!(ys,maximum.(cumsum.(read_inf_data(g,dloc=dloc,beta=beta,gamma=gamma,dtype="cinfs"))./nnodes))
#     end
#
#
#     f1 = scatter(1:length(ys[1]),ys[1],ylabel="Peak Infections",ylims=(0,1),
#         xticks=(250:500:8000,string.(0:16)),
#         xlabel="quarantine percent",label="",
#         title="$(gnames[1]) - Maximum Infections\n$(uppercase(method))($beta,$gamma)")
#
#     f2 = scatter(1:length(ys[2]),ys[2],ylabel="Peak Infections",ylims=(0,1),
#         xticks=(250:500:8000,string.(0:16)),
#         xlabel="quarantine percent",label="",
#         title="$(gnames[2]) - Maximum Infections\n$(uppercase(method))($beta,$gamma)")
#
#     f3 = scatter(1:length(ys[3]),ys[3],ylabel="Peak Infections",ylims=(0,1),
#         xticks=(250:500:8000,string.(0:16)),
#         xlabel="quarantine percent",label="",
#         title="$(gnames[3]) - Maximum Infections\n$(uppercase(method))($beta,$gamma)")
#
#     Plots.savefig(f1,dst*gnames[1]*"-$beta-max-infs.png")
#     Plots.savefig(f2,dst*gnames[2]*"-$beta-max-infs.png")
#     Plots.savefig(f3,dst*gnames[3]*"-$beta-max-infs.png")
#
#     return f1,f2,f3
# end


# for h in gs
#     gname = getgnames(h,gpath)[1][1:end-5]
#     nnodes = size(loadGraph(gname*".smat",gpath),1)
#     for beta in [collect(0.001:0.001:0.01); collect(0.02:0.01:0.1)]
#         println("working on beta $beta for $gname")
#         peakInfsCompPlot(gname,7,dloc=dloc,dst=dst,beta=beta,nnodes=nnodes)
#     end
# end
#

# peakInfsCompPlot(getgnames(gs[1],gpath)[1][1:end-5],7;gpath=gpath,dloc=dloc,dst="purgatory/")

## sparsification heatmaps

# function graph_sparification_heatmap(gname::String,beta::Float64,gamma::Float64=0.05;
#                         gpath::String="",dtype::String="tinfs",dloc::String="",
#                         method::String="sir",
#                         rps::Vector{Float64}=collect(0.1:0.1:0.9))
#
#     gr()
#     # do stuff
#     A = loadGraph(gname,gpath)
#     nnodes = size(A,1)
#     ths = nnodes.*collect(0.01:0.01:1.0)
#     g = gname[1:end-5]
#
#     data = Vector{Vector{Vector{Int}}}()
#     gnames = Vector{String}()
#     push!(data,cumsum.(read_inf_data(g,dloc=dloc,beta=beta,gamma=gamma,dtype=dtype)))
#     push!(gnames,g)
#     for p in reverse(rps)
#         push!(data,cumsum.(read_inf_data("cn-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype)))
#         push!(gnames,"cn-$p-rewired-$g")
#     end
#
#     #build heatmap
#     #will do color as pr( >=k% of nodes)
#     pdata = zeros(length(ths),length(gnames))
#     for graphind=1:length(gnames)
#         if dtype == "tinfs"
#             ps = map(x->findlast(x.>ths),last.(data[graphind][1:500]))
#         elseif dtype == "cinfs"
#             ps = map(x->findlast(x.>ths),maximum.(data[graphind][1:500]))
#         end
#         ps = Vector{Union{Int,Nothing}}(ps)
#         ps[ps.==nothing].= 0
#         ps = Int.(ps)
#         pdata[:,graphind] .= fit(Histogram,ps,1:length(ths)+1).weights
#         pdata[:,graphind] .= cumsum(pdata[end:-1:1,graphind])[end:-1:1]./500 #500 was the sample size for the data
#     end
#     if dtype == "tinfs"
#         ylab = "Pr(Total Infs ≥ k% of Nodes)"
#     elseif dtype == "cinfs"
#         ylab = "Pr(Maximum Infs ≥ k% of Nodes)"
#     end
#     f1 = heatmap(pdata,ylabel="Percent of Total Nodes",colorbar_title=ylab,
#                 title="$(gname[1:end-5]) - SIR($beta,$gamma)",clims=(0,1),
#                 xlabel="     Original ⇒ Sparsified Graph\nSparsification Parameter",
#                 xrotation = 90)
#
#     if dtype == "tinfs"
#         clabel = "Fraction of Infected Nodes"
#         data = hcat(map(x->last.(x),data)...)./nnodes
#     elseif dtype == "cinfs"
#         clabel = "Normalized Maximum Infection Size"
#         data = hcat(map(x->maximum.(x),data)...)./nnodes
#     end
#     f = heatmap(data,ylabel="quarantine percentage",
#         xlabel="     Original ⇒ Sparsified Graph\nSparsification Parameter",
#         title="$(gname[1:end-5]) - SIR($beta,$gamma)\n500 node sample",
#         colorbar_title=clabel,c=:heat,clims=(0,1),
#         yticks=(collect(250:500:8000),0:15),xrotation = 90)
#
#     plot!(f,xticks=(1:length(rps)+1, string.(round.(100 .*vec([1.0 reverse(rps)...]),digits=2))),margins=9mm)
#     plot!(f1,xticks=(1:length(rps)+1, string.(round.(100 .*vec([1.0 reverse(rps)...]),digits=2))),margins=9mm)
#     f,f1
# end


# gpath = "graphs/all-graphs/"
# dloc = "diffusion-data/tmp/"
# ps = collect(0.1:0.1:0.9)
#
# gs = ["astro","enron","ucf","unc","penn"]
# bs = sort([collect(0.001:0.001:0.01);collect(0.01:0.01:0.1)])#bs])
#
# gnames = getgnames(h,gpath)
# filter!(x->!startswith(x,"er-") && !startswith(x,"rewired-") && !startswith(x,"cl-"),gnames)
# for gname in [gnames[1]]
#     println("working on $gname")
#     for beta in bs
#         #total infs plots
#         f,f1 = graph_sparification_heatmap(gname,beta,0.05,gpath=gpath,dloc=dloc,dtype="tinfs",rps = ps)


# for h in gs
#     gnames = getgnames(h,gpath)
#     filter!(x->!startswith(x,"er-") && !startswith(x,"rewired-") && !startswith(x,"cl-"),gnames)
#     for gname in [gnames[1]]
#         println("working on $gname")
#         for beta in bs
#             #total infs plots
#             f,f1 = graph_sparification_heatmap(gname,beta,0.05,gpath=gpath,dloc=dloc,dtype="tinfs",rps = ps)
#             Plots.savefig(f,"purgatory/totalinfs-sparsification-$(gname[1:end-5])-$beta-1.png")
#             Plots.savefig(f1,"purgatory/totalinfs-sparsification-$(gname[1:end-5])-$beta-2.png")
#             #maximum infs plots
#             f,f1 = graph_sparification_heatmap(gname,beta,0.05,gpath=gpath,dloc=dloc,dtype="cinfs",rps = ps)
#             Plots.savefig(f,"purgatory/peakinfs-sparsification-$(gname[1:end-5])-$beta-1.png")
#             Plots.savefig(f1,"purgatory/peakinfs-sparsification-$(gname[1:end-5])-$beta-2.png")
#         end
#     end
# end


# g = getgnames("astro","graphs/all-graphs/")[1]
# f,f1 = graph_sparification_heatmap(g,0.001,gpath="graphs/all-graphs/",dloc="diffusion-data/tmp/")
#
# f
# f1




# include("fast-diffusion.jl")
#
# using Measures
# function graph_heatmap1(gname::String,beta::Float64,gamma::Float64=0.05;
#                         gpath::String="",dtype::String="tinfs",dloc::String="",
#                         dloc2::String="diffusion-data/tmp/",method::String="sir",
#                         rps::Vector{Float64}=vec([0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.25 0.5 1.0 10.0]))
#
#     gr()
#     # do stuff
#     A = loadGraph(gname,gpath)
#     nnodes = size(A,1)
#     ths = nnodes.*collect(0.01:0.01:1.0)
#     g = gname[1:end-5]
#
#     data = Vector{Vector{Vector{Int}}}()
#     gnames = Vector{String}()
#     for p in 100 .*reverse(rps)# 0.1 0.25 0.5 1.0]# 10.0]#[1 5 10 25 50 100 1000]
#         a = [cumsum.(read_inf_data("rewired-$p-$g",dloc=dloc2,beta=beta,gamma=gamma,dtype=dtype))[1:500];cumsum.(read_inf_data("rewired-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype))]
#         push!(data,a)
#         push!(gnames,"rewired-$p-$g")
#     end
#     push!(data,[cumsum.(read_inf_data(g,dloc=dloc2,beta=beta,gamma=gamma,dtype=dtype))[1:500];cumsum.(read_inf_data(g,dloc=dloc,beta=beta,gamma=gamma,dtype=dtype))])
#     # push!(data,cumsum.(read_inf_data(g,dloc=dloc,beta=beta,gamma=gamma,dtype=dtype)))
#     push!(gnames,g)
#     for p in 100 .*rps# 0.1 0.25 0.5 1.0]# 10.0]#[1 5 10 25 50 100 1000]
#         push!(data,[cumsum.(read_inf_data("er-$p-$g",dloc=dloc2,beta=beta,gamma=gamma,dtype=dtype))[1:500];cumsum.(read_inf_data("er-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype))])
#         # push!(data,cumsum.(read_inf_data("er-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype)))
#         # push!(data,cumsum.(read_inf_data("er-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype)))
#         push!(gnames,"er-$p-rewired-$g")
#     end
#
#     #build heatmap
#     #will do color as pr( >=k% of nodes)
#     pdata = zeros(length(ths),length(gnames))
#     for graphind=1:length(gnames)
#         if dtype == "tinfs"
#             ps = map(x->findlast(x.>ths),last.(data[graphind][1:500]))
#         elseif dtype == "cinfs"
#             ps = map(x->findlast(x.>ths),maximum.(data[graphind][1:500]))
#         end
#         ps = Vector{Union{Int,Nothing}}(ps)
#         ps[ps.==nothing].= 0
#         ps = Int.(ps)
#         pdata[:,graphind] .= fit(Histogram,ps,1:length(ths)+1).weights
#         pdata[:,graphind] .= cumsum(pdata[end:-1:1,graphind])[end:-1:1]./500 #500 was the sample size for the data
#     end
#     if dtype == "tinfs"
#         ylab = "Pr(Total Infs ≥ k% of Nodes)"
#     elseif dtype == "cinfs"
#         ylab = "Pr(Maximum Infs ≥ k% of Nodes)"
#     end
#     f1 = heatmap(pdata,ylabel="Percent of Total Nodes",colorbar_title=ylab,
#                 title="$(gname[1:end-5]) - SIR($beta,$gamma)",clims=(0,1),
#                 xlabel="     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent",
#                 xrotation = 90)
#
#     if dtype == "tinfs"
#         clabel = "Fraction of Infected Nodes"
#         data = hcat(map(x->last.(x),data)...)./nnodes
#     elseif dtype == "cinfs"
#         clabel = "Normalized Maximum Infection Size"
#         data = hcat(map(x->maximum.(x),data)...)./nnodes
#     end
#     f = heatmap(data,ylabel="rolling quarantine",
#         xlabel="     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent",
#         title="$(gname[1:end-5]) - SIR($beta,$gamma)\n500 node sample",
#         colorbar_title=clabel,c=:heat,clims=(0,1),
#         yticks=(collect(250:500:5500),0:100:1000),xrotation = 90)
#
#     plot!(f,xticks=(1:2*length(rps)+1, string.(round.(100 .*vec([reverse(rps)... 0 rps...]),digits=2))),margins=9Measures.mm)
#     plot!(f1,xticks=(1:2*length(rps)+1, string.(round.(100 .*vec([reverse(rps)... 0 rps...]),digits=2))),margins=9Measures.mm)
#
#     return f,f1
# end

# ps = vec([0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.25 0.5 1.0 10.0 100.0])
#
# gs = ["cn-modBerk","cn-modStan","cn-penn"]
# gs = ["anony"]
# for gname in gs
#     g = getgnames(gname,"graphs/all-graphs/")[1]
#     println("working on $g")
#     for b in [collect(0.001:0.001:0.01);collect(0.02:0.01:0.1)]
#         println("working on $b")
#         f,f1 = graph_heatmap1(g,b,gpath="graphs/all-graphs/",dloc="diffusion-data/qtmp/")
#         Plots.savefig(f,"purgatory/totalinfs-fixed-q-$(g[1:end-5])-$b-1.png")
#         Plots.savefig(f1,"purgatory/totalinfs-fixed-q-$(g[1:end-5])-$b-2.png")
#     end
# end


# getgnames("anony","graphs/all-graphs/")[1]

# function graph_heatmap_high_percentage(gname::String,beta::Float64,gamma::Float64=0.05;
#                         gpath::String="",dtype::String="tinfs",dloc::String="",
#                         dloc2::String="diffusion-data/tmp/",method::String="sir",
#                         rps::Vector{Float64}=vec([0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.25 0.5 1.0 10.0]))
#
#     gr()
#     # do stuff
#     A = loadGraph(gname,gpath)
#     nnodes = size(A,1)
#     ths = nnodes.*collect(0.01:0.01:1.0)
#     g = gname[1:end-5]
#
#     data = Vector{Vector{Vector{Int}}}()
#     gnames = Vector{String}()
#     for p in 100 .*reverse(rps)# 0.1 0.25 0.5 1.0]# 10.0]#[1 5 10 25 50 100 1000]
#         a = [cumsum.(read_inf_data("rewired-$p-$g",dloc=dloc2,beta=beta,gamma=gamma,dtype=dtype))[vcat(1:500,2501:3000,5001:5500,7501:8000)];cumsum.(read_inf_data("rewired-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype))]
#         push!(data,a)
#         push!(gnames,"rewired-$p-$g")
#     end
#     push!(data,[cumsum.(read_inf_data(g,dloc=dloc2,beta=beta,gamma=gamma,dtype=dtype))[vcat(1:500,2501:3000,5001:5500,7501:8000)];cumsum.(read_inf_data(g,dloc=dloc,beta=beta,gamma=gamma,dtype=dtype))])
#     # push!(data,cumsum.(read_inf_data(g,dloc=dloc,beta=beta,gamma=gamma,dtype=dtype)))
#     push!(gnames,g)
#     for p in 100 .*rps# 0.1 0.25 0.5 1.0]# 10.0]#[1 5 10 25 50 100 1000]
#         push!(data,[cumsum.(read_inf_data("er-$p-$g",dloc=dloc2,beta=beta,gamma=gamma,dtype=dtype))[vcat(1:500,2501:3000,5001:5500,7501:8000)];cumsum.(read_inf_data("er-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype))])
#         # push!(data,cumsum.(read_inf_data("er-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype)))
#         # push!(data,cumsum.(read_inf_data("er-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype)))
#         push!(gnames,"er-$p-rewired-$g")
#     end
#
#     #build heatmap
#     #will do color as pr( >=k% of nodes)
#     # pdata = zeros(length(ths),length(gnames))
#     # for graphind=1:length(gnames)
#     #     if dtype == "tinfs"
#     #         ps = map(x->findlast(x.>ths),last.(data[graphind][1:500]))
#     #     elseif dtype == "cinfs"
#     #         ps = map(x->findlast(x.>ths),maximum.(data[graphind][1:500]))
#     #     end
#     #     ps = Vector{Union{Int,Nothing}}(ps)
#     #     ps[ps.==nothing].= 0
#     #     ps = Int.(ps)
#     #     pdata[:,graphind] .= fit(Histogram,ps,1:length(ths)+1).weights
#     #     pdata[:,graphind] .= cumsum(pdata[end:-1:1,graphind])[end:-1:1]./500 #500 was the sample size for the data
#     # end
#     if dtype == "tinfs"
#         ylab = "Pr(Total Infs ≥ k% of Nodes)"
#     elseif dtype == "cinfs"
#         ylab = "Pr(Maximum Infs ≥ k% of Nodes)"
#     end
#     # f1 = heatmap(pdata,ylabel="Percent of Total Nodes",colorbar_title=ylab,
#     #             title="$(gname[1:end-5]) - SIR($beta,$gamma)",clims=(0,1),
#     #             xlabel="     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent",
#     #             xrotation = 90)
#
#     if dtype == "tinfs"
#         clabel = "Fraction of Infected Nodes"
#         data = hcat(map(x->last.(x),data)...)./nnodes
#     elseif dtype == "cinfs"
#         clabel = "Normalized Maximum Infection Size"
#         data = hcat(map(x->maximum.(x),data)...)./nnodes
#     end
#     f = heatmap(data,ylabel="rolling quarantine",
#         xlabel="     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent",
#         title="$(gname[1:end-5]) - SIR($beta,$gamma)\n500 node sample",
#         colorbar_title=clabel,c=:heat,clims=(0,1),
#         yticks=(collect(250:500:8000),0:5:75),xrotation = 90)
#
#     plot!(f,xticks=(1:2*length(rps)+1, string.(round.(100 .*vec([reverse(rps)... 0 rps...]),digits=2))),margins=9Measures.mm)
#     # plot!(f1,xticks=(1:2*length(rps)+1, string.(round.(100 .*vec([reverse(rps)... 0 rps...]),digits=2))),margins=9Measures.mm)
#
#     return f
# end

# dst = "purgatory/"
# gs = ["mexico","dblp","anony"]
# for g in gs
#     println("working on $g")
#     g = getgnames(g,"graphs/all-graphs/")[1]
#     for b in [collect(0.001:0.001:0.01);collect(0.01:0.01:0.1)]
#         f = graph_heatmap_high_percentage(g,b,gpath="graphs/all-graphs/",dloc="diffusion-data/high-percentage/")
#         Plots.savefig(f,dst*"totalinfs-high-percentage-$(g[1:end-5])-$b.png")
#     end
# end






# using Plots
# using Measures
# function single_graph_heatmap(gname::String;
#                                 betas::Vector{Float64}=[collect(0.001:0.001:0.01);collect(0.02:0.01:0.1)],
#                                 gamma::Float64=0.05,
#                                 gpath::String="",
#                                 dtype::String="tinfs",
#                                 dloc::String="",
#                                 method::String="sir")
#
#     #start of fcn
#     gr()
#     A = loadGraph(gname,gpath)
#     nnodes = size(A,1)
#     ths = nnodes.*collect(0.01:0.01:1.0)
#     g = gname[1:end-5]
#
#     data = Vector{Vector{Vector{Int}}}()
#     for b in betas
#         push!(data,cumsum.(read_inf_data("$g",dloc=dloc,beta=b,gamma=gamma,dtype=dtype)))
#     end
#
#     #build heatmap
#     #will do color as pr( >=k% of nodes)
#     pdata = zeros(length(ths),length(betas))
#     for bind=1:length(betas)
#         if dtype == "tinfs"
#             ps = map(x->findlast(x.>ths),last.(data[bind][1:500]))
#         elseif dtype == "cinfs"
#             ps = map(x->findlast(x.>ths),maximum.(data[bind][1:500]))
#         end
#         ps = Vector{Union{Int,Nothing}}(ps)
#         ps[ps.==nothing].= 0
#         ps = Int.(ps)
#         pdata[:,bind] .= fit(Histogram,ps,1:length(ths)+1).weights
#         pdata[:,bind] .= cumsum(pdata[end:-1:1,bind])[end:-1:1]./500 #500 was the sample size for the data
#     end
#     if dtype == "tinfs"
#         ylab = "Pr(Total Infs ≥ k% of Nodes)"
#     elseif dtype == "cinfs"
#         ylab = "Pr(Maximum Infs ≥ k% of Nodes)"
#     end
#     f1 = heatmap(pdata,ylabel="Percent of Total Nodes",colorbar_title=ylab,
#                 title="$(gname[1:end-5]) - SIR(β,$gamma)",clims=(0,1),
#                 xlabel="β",xrotation = 90)
#
#     if dtype == "tinfs"
#         clabel = "Fraction of Infected Nodes"
#         data = hcat(map(x->last.(x),data)...)./nnodes
#     elseif dtype == "cinfs"
#         clabel = "Normalized Maximum Infection Size"
#         data = hcat(map(x->maximum.(x),data)...)./nnodes
#     end
#     f = heatmap(data,ylabel="quarantine percentage",
#         xlabel="β",
#         title="$(g) - SIR(β,$gamma)\n500 node sample",
#         colorbar_title=clabel,c=:heat,clims=(0,1),
#         yticks=(collect(250:500:8000),0:15),xrotation = 90)
#
#     plot!(f,xticks=(1:length(betas), string.(betas)),margins=9Measures.mm)
#     plot!(f1,xticks=(1:length(betas), string.(betas)),margins=9Measures.mm)
#
#     return f,f1
# end

# gpath = "graphs/all-graphs/"
# gpath = "graphs/synthetic-graphs/"
# gnames = getgnames("rw",gpath)
# f,f1 = single_graph_heatmap(gnames[1],gpath=gpath,dloc="../../qtmp/")

# gpath = "purgatory/diff/diffusion-graphs/"
# dst = "purgatory/heatmaps/"
# for g in gnames[2:end]
# gnames = getgnames("lfr",gpath)
# g = gnames[4]
# make_ncp_plots(g,gpath=gpath,dst=dst)
#
#     f,f1 = single_graph_heatmap(g,gpath=gpath,dloc=dloc)
#
#     Plots.savefig(f,dst*"tinfs-beta-heatmap-q-$(g[1:end-5])-sir-0.05.png")
#     Plots.savefig(f1,dst*"tinfs-beta-heatmap-$(g[1:end-5])-sir-0.05.png")
# end
# f
# f1
#
# #reproducing lldm figures
# gpath = "graphs/synthetic-graphs/lldm-forest-fire/"
# dst = "purgatory/ncps/forest-fire/"
# gnames = getgnames("forest",gpath)
# for gname in gnames
#     A = loadGraph(gname,gpath)
#     # A = largest_component(A)[1]
#     davg = nnz(A)/size(A,1)
#     ncp,sets = make_ncp(A,get_sets=true)
#     f = plot_ncp(ncp,A)
#     title!(f,"NCP - $(gname[1:end-5])\ndavg: $(round(davg,digits=2)), nnodes: $(size(A,1))")
#     xlims!(f,(1,1e4))
#     ylims!(f,(1e-3,1))
#     Plots.savefig(f,dst*"ncp-$(gname[1:end-5]).png")
#
#     f = ncpcoverage(ncp,sets,A)
#     title!(f,"Coverage - $(gname[1:end-5])\ndavg: $(round(davg,digits=2)), nnodes: $(size(A,1))")
#     xlims!(f,(1,1e4))
#     ylims!(f,(1e-3,1))
#     Plots.savefig(f,dst*"ncpcoverage-$(gname[1:end-5]).png")
# end
#
#
#
#
# # gname = gnames[55]
# # dloc ="diffusion-data/tmp/"
#
# dloc = "diffusion-data/synthetic-graphs/"
# for (i,g) in enumerate(gnames)
#     println("$i of $(length(gnames))")
#     f,f1 = single_graph_heatmap(g,gpath=gpath,dloc=dloc)
#     #save images
#
#     Plots.savefig(f,"purgatory/heatmaps/tinfs-beta-heatmap-q-$(g[1:end-5])-sir-0.05.png")
#     Plots.savefig(f1,"purgatory/heatmaps/tinfs-beta-heatmap-$(g[1:end-5])-sir-0.05.png")
# end


## scatter plots
# include("fast-diffusion.jl")
# include("ncp-computation-script.jl")
# gnames = getgnames("lfr","graphs/synthetic-graphs/lfr-graphs/")
# g = gnames[70]
# f,f1 = single_graph_heatmap(g,gpath="graphs/synthetic-graphs/lfr-graphs/",dloc="diffusion-data/synthetic-graphs/")
# f1
# f
# gnames = getgnames("mexico","graphs/all-graphs/")
# g = gnames[1]
# dloc = "diffusion-data/tmp/"
# a = read_inf_data(g[1:end-5],dloc=dloc,beta=0.05,dtype="tinfs")
#
# ys = cumsum.(a[1:500])
# inds = map(x->findfirst(x.>=10000),ys)
# inds = inds[deleteat!(collect(1:length(inds)),inds.==nothing)]
# histogram(inds,xlims=(1,200),xlabel="time to 10k infs",ylabel="count",leg=false)
#
#
# gpath = "graphs/synthetic-graphs/lfr-graphs/"
# gnames = getgnames("lfr",gpath)
#
# gname = gnames[64]
# A = MatrixNetworks.readSMAT(gpath*"lfr-100000-2.30-1.50-0.01.smat")
#
# A = max.(A,A')
# nnz(A)
#
#
# hnames = readdir("graphs/all-graphs/")
# gnames = readdir("graphs/rewired-graphs/")
# ncp stuff
# dloc = "ncpdata/"
# gpath = "graphs/all-graphs/"
# gs = ["anony","dblp","mexico"]
# gs = ["study","wisconsin","commutes",
#         "usf","ucf","uf","astro"]
# for g in gs
#     gnames = getgnames(g,gpath)
#     filter!(x->!startswith(x,"cl-") && !startswith(x,"rewired-") && !startswith(x,"er-"),gnames)
#     #load ncp sests
#     for gname in gnames
#         A = loadGraph(gname,gpath)
#         ncp,headerinfo = readdlm(dloc*"ncpinfo-$(gname[1:end-5]).txt",',',header=true)
#         ncp = DataFrame(ncp,vec(headerinfo))
#         ncpmat = MatrixNetworks.readSMAT(dloc*"ncpmat-$gname")
#         sets = ncpmat_to_sets(ncpmat)
#
#         make_ncp_plots(ncp,sets,A,gname=gname,dst="purgatory/")
#     end
# end
# using DataFrames
#
#
#
# gnames = getgnames("mex",gpath)
# filter!(x->!startswith(x,"cl-") && !startswith(x,"rewired-") && !startswith(x,"er-"),gnames)
# gname = gnames[1]
# #load ncp sests
# A = loadGraph(gname,gpath)
# ncp,headerinfo = readdlm(dloc*"ncpinfo-$(gname[1:end-5]).txt",',',header=true)
# ncp = DataFrame(ncp,vec(headerinfo))
# ncpmat = MatrixNetworks.readSMAT(dloc*"ncpmat-$gname")
# sets = ncpmat_to_sets(ncpmat)
# pyplot()
#
# ncp.size = map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1))
# myncpplot(ncp.size,ncp.cond)
#
# x,y = ncp.size,ncp.cond
# nbins = 100
# myhexbin(x,y,nbins=nbins)
#
#
#
# hexhist = fit(HexBinPlots.HexHistogram,log10.(x),log10.(y),nbins)
# h,vh = HexBinPlots.make_shapes(hexhist)
# vmax = maximum(vh)
# xsh = Vector{Float64}()
# ysh = Vector{Float64}()
# for k in eachindex(h)
#     append!(xsh,h[k].x)
#     push!(xsh,h[k].x[1])
#     push!(xsh,NaN)
#     append!(ysh,h[k].y)
#     push!(ysh,h[k].y[1])
#     push!(ysh,NaN)
# end
# # @show vh #no need to display this
# pyplot()
# Plots.plot(10.0.^xsh,10.0.^ysh,fill_z=log10.(vh.+1),linecolor=nothing,
#     seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=false)


##
#lfr ncps w/ planted communities
#planted communities labels not the same since we look at largest component
# using MatrixNetworks
# using DelimitedFiles
# using SparseArrays
# gpath = "graphs/synthetic-graphs/lfr-graphs/new-graphs/"
# ngpath = "graphs/tmpgraphs1/"
#
# g = "lfr-100000-3.00-2.00-0.15-1-2000-5-500-17"
# cmap = vec(readdlm(gpath*"cmap-$(g).txt",Int))
#
# A = MatrixNetworks.readSMAT(gpath*g*".smat")
# A = max.(A,A')
# fill!(A.nzval,1.0)
# da = vec(sum(A;dims=2))
# function conductance(A::SparseMatrixCSC,S::Set{Int};rowids = rowvals(A))
#     vol = 0.0
#     bd = 0.0
#     @inbounds for v in S
#         for u in rowids[nzrange(A,v)]
#             if !(u in S) #cut edge
#                 bd += 1
#             end
#             vol += 1
#         end
#     end
#     vol = min(vol,nnz(A)-vol)
#     return bd/vol
# end
# function build_sets(cmap::Vector{Int})
#     sets = Dict{Int,Set{Int}}()
#     for (u,c) in enumerate(cmap)
#         if haskey(sets,c)
#             push!(sets[c],u)
#         else
#             sets[c] = Set(u)
#         end
#     end
#     return sets
# end
#
# #when we connected graphs, we sorted nodes by degree so need to undo that
# wvec = vec(sum(A;dims=2))
# nodeperm = sortperm(wvec,rev=true)
# pinverse = sortperm(nodeperm)
#
# gnames = getgnames(g,ngpath)
# filter!(x->occursin("connect-graph",x),gnames)
# for gname in gnames
#     #make ncp plots
#     println("working on $gname")
#     X = MatrixNetworks.readSMAT(ngpath*gname)
#     X = max.(X,X')-Diagonal(X)
#     fill!(X.nzval,1.0)
#     # X = X[pinv,pinv]
#     X,bmap = largest_component(X)
#
#     newcmap = cmap[nodeperm][bmap]
#
#     open(ngpath*"cmap-gc-$(gname[1:end-5]).txt","w") do io
#         writedlm(io,newcmap)
#     end
# end
#=
    szs = Dict(k=>length(sets[k]) for k in keys(sets))

    rowids = rowvals(X)
    conds = map(x->conductance(X,sets[x],rowids=rowids),collect(keys(sets)))
    szs = collect(values(szs))

    #do ncp stuff
    ncp,headerinfo = readdlm(dloc*"ncpinfo-$(gname[1:end-5]).txt",',',header=true)
    ncp = DataFrame(ncp,vec(headerinfo))
    ncpmat = MatrixNetworks.readSMAT(dloc*"ncpmat-$gname")
    sets = ncpmat_to_sets(ncpmat)

    f = make_hexbinncp_plots(ncp,sets,X,gname=g,dst="purgatory/ncps/")
    scatter!(f,szs,conds,marker=:star,marker_size=1,c=:red,label="planted\ncommunities")
    Plots.savefig(f,dst*"hexbinncp-planted-$(gname[1:end-5]).png")
end
=#

#do original graph
# X,bmap = largest_component(A)
# newcmap = cmap[bmap]
# open(gpath*"cmap-gc-$g.txt","w") do io
#     writedlm(io,newcmap)
# end


##

#=
#method 1
# undo node perm then find new labels in largest component of lfr graph
gname = gnames[1]
fname = "cmap-gc-"*gname[1:end-5]*".txt"
ncmap = vec(readdlm(ngpath*fname,Int))
#make comms
sets = Dict{Int,Set{Int}}()
for (u,c) in enumerate(ncmap)
    if haskey(sets,c)
        push!(sets[c],u)
    else
        sets[c] = Set(u)
    end
end

X = MatrixNetworks.readSMAT(ngpath*gname)
X = max.(X,X')-Diagonal(X)
fill!(X.nzval,1.0)
X = X[pinv,pinv]
X,bmap = largest_component(X)

rowids = rowvals(X)
map(x->conductance(X,sets[x],rowids=rowids),collect(keys(sets)))
=#

#method 2
#apply node perm then take largest component #TODO this is the one we want!!!
# gname = gnames[1]
# fname = "cmap-gc-"*gname[1:end-5]*".txt"
# ncmap = vec(readdlm(ngpath*fname,Int))
#
# A = MatrixNetworks.readSMAT(ngpath*gname)
# A = max.(A,A')-Diagonal(A)
# fill!(A.nzval,1.0)
# A,bmap = largest_component(A)
#
# newcmap = cmap[nodeperm][bmap]
# sets = Dict{Int,Set{Int}}()
# for (u,c) in enumerate(newcmap)
#     if haskey(sets,c)
#         push!(sets[c],u)
#     else
#         sets[c] = Set(u)
#     end
# end
#
# rowids = rowvals(A)
# map(x->conductance(A,sets[x],rowids=rowids),collect(keys(sets)))

#NEW STUFF
#original graph

# gpath = "graphs/synthetic-graphs/lfr-graphs/new-graphs/"
# ngpath = "graphs/tmpgraphs1/"
#
# g = "lfr-100000-3.00-2.00-0.15-1-2000-5-500-17"
# cmap = vec(readdlm(gpath*"cmap-$(g).txt",Int))
#
# A = MatrixNetworks.readSMAT(gpath*g*".smat")
# A = max.(A,A')
# fill!(A.nzval,1.0)
# A = A[nodeperm,nodeperm]
# Random.seed!(0)
# A = fast_chung_lu_connect_components(A)
# A = A[pinverse,pinverse]
# A,bmap = largest_component(A)
#
# sets = build_sets(cmap[bmap])
# rowids = rowvals(A)
# map(x->conductance(A,sets[x],rowids=rowids),collect(keys(sets)))
#
#
#
#
# ## lfr planted comms comparison
# include("fast-diffusion.jl")
# using Plots
# using Interact
# using MatrixNetworks
# using DelimitedFiles
# using SparseArrays
# gpath = "graphs/lfr-graphs/"
#
#
# function build_sets(cmap::Vector{Vector{Int}})
#     sets = Dict{Int,Set{Int}}()
#     for (node,cs) in enumerate(cmap)
#         for c in cs
#             if haskey(sets,c)
#                 push!(sets[c],node)
#             else
#                 sets[c] = Set(node)
#             end
#         end
#     end
#     return sets
# end
# function read_cmap(fname::String)
#     result = Vector{Vector{Int}}()
#     open(fname,"r") do io
#         while !eof(io)
#             newline = tryparse.(Int,split(readline(io),","))
#             if newline[1] == nothing
#                 push!(result,[-1])
#             else
#                 push!(result,newline)
#             end
#         end
#     end
#     return result
# end
# function conductance(A::SparseMatrixCSC,S::Set{Int};rowids = rowvals(A))
#     vol = 0.0
#     bd = 0.0
#     @inbounds for v in S
#         for u in rowids[nzrange(A,v)]
#             if !(u in S) #cut edge
#                 bd += 1
#             end
#             vol += 1
#         end
#     end
#     vol = min(vol,nnz(A)-vol)
#     return bd/vol
# end
#
# g = "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17"
# gnames = getgnames(g,gpath)
# nsamples = [0;100; collect(500:500:8000)]
# # nsamples = collect(500:500:8000)
# #load graph
# #load maps
# #convert maps to sets
# #get conductances
# #compare planted sets to sets after refining
#
# #TODO just missing a couple of files here.
# @manipulate throttle=1 for numRWs in nsamples
#     if numRWs==0
#         gname = g*".smat"
#     else
#         gname = g*"-connect-graph-invdegdepth-"*string(numRWs)*"-0.5-0.0-100-5.smat"
#     end
#     A = loadGraph(gname,gpath)
#     cmap = read_cmap(gpath*"cmap-gc-$(gname[1:end-5]).txt")
#     mqimap = read_cmap(gpath*"cmap-gc-mqi-$(gname[1:end-5]).txt")
#     slmap = read_cmap(gpath*"cmap-gc-sl-$(gname[1:end-5]).txt")
#
#     sets = build_sets(cmap)
#     mqisets = build_sets(mqimap)
#     slsets = build_sets(slmap)
#
#     ks = sort(collect(keys(sets)))
#
#     rowids = rowvals(A)
#     c1 = map(x->conductance(A,sets[x],rowids=rowids),ks)
#     c2 = map(x->conductance(A,mqisets[x],rowids=rowids),ks)
#     c3 = map(x->conductance(A,slsets[x],rowids=rowids),ks)
#     # inds = findall(all.(zip(abs.(c1.-c3).>=0.3,c3.<0.1)))
#     inds = findall(c3.<=0.05)
#     ls = map(x->length(sets[x]),ks[inds])
#     # delta = map(x->length(symdiff(sets[x],slsets[x])),ks[inds])
#     delta = map(x->length(intersect(sets[x],slsets[x]))/length(sets[x]),ks[inds])
#
#     p1 = scatter(inds,c1[inds],yscale=:log10,ylims=(1e-6,1.1),label="planted set",legend=:bottomleft,c=:black,xlims=(0,length(ks)+1))
#     # scatter!(p1,inds,c2[inds],label="MQI")
#     scatter!(p1,inds,c3[inds],label="SimpleLocal")
#     plot!(p1,ylabel="Conductance",
#             xlabel="Set",
#             title="LFR+RW - RW samples: $numRWs\nFinal Conductance ≤ 0.05\nTotal Sets: $(length(inds))")
#     # histogram!(ls,bins=0:10:100,ylabel="",xlabel="size of planted set",label="",title="Total Sets: $(length(ls))",
#     #     xlims=(0,100),ylims=(0,500),
#     #     titlefontsize=10,
#     #     inset=bbox(0.7, 0.21, 0.2, 0.2, :bottom),
#     #     subplot=2)
#     # histogram!(delta,bins=0:10:100,ylabel="",xlabel="size of sym. difference",label="",title="Total Sets: $(length(delta))",
#     #     xlims=(0,100),ylims=(0,500),
#     #     titlefontsize=10,
#     #     inset=bbox(0.7, 0.215, 0.2, 0.2, :bottom),
#     #     subplot=2)
#     histogram!(ls,bins=0:10:100,
#         normalize=:probability,
#         ylabel="",
#         xlabel="|Planted|",
#         label="",
#         xlims=(0,100),ylims=(-0.01,0.4),
#         xguidefontsize=10,
#         inset=bbox(0.4, 0.21, 0.2, 0.2, :bottom),
#         subplot=2)
#     histogram!(delta,bins=0:0.1:1.0,ylabel="",
#         xlabel="|Planted ⋂ New|/|Planted|",
#         label="",
#         normalize=:probability,
#         # title="Planted vs SimpleLocal Sets",
#         xlims=(0,1.0),ylims=(-0.01,0.4),
#         xguidefontsize=10,
#         inset=bbox(0.72, 0.215, 0.2, 0.2, :bottom),
#         subplot=3,
#         c=:green)
# end



# for nsample in nsamples
#     gname = g*"-connect-graph-invdegdepth-"*string(nsample)*"-0.5-0.0-100-5.smat"
#     A = loadGraph(gname,gpath)
#     cmap = read_cmap(gpath*"cmap-gc-$(gname[1:end-5]).txt")
#     mqimap = read_cmap(gpath*"cmap-gc-mqi-$(gname[1:end-5]).txt")
#     slmap = read_cmap(gpath*"cmap-gc-sl-$(gname[1:end-5]).txt")
#
#     sets = build_sets(cmap)
#     mqisets = build_sets(mqimap)
#     slsets = build_sets(slmap)
#
#     ks = sort(collect(keys(sets)))
#
#     rowids = rowvals(A)
#     c1 = map(x->conductance(A,sets[x],rowids=rowids),ks)
#     c2 = map(x->conductance(A,mqisets[x],rowids=rowids),ks)
#     c3 = map(x->conductance(A,slsets[x],rowids=rowids),ks)
#     # inds = findall(all.(zip(abs.(c1.-c3).>=0.3,c3.<0.1)))
#     inds = findall(c3.<0.05)
#     ls = map(x->length(sets[x]),ks[inds])
#
#     p1 = scatter(inds,c1[inds],yscale=:log10,ylims=(1e-6,1.1),label="planted set",legend=:bottomleft,c=:black,xlims=(0,length(ks)+1))
#     # scatter!(p1,inds,c2[inds],label="MQI")
#     scatter!(p1,inds,c3[inds],label="SimpleLocal")
#     plot!(p1,ylabel="conductance",xlabel="planted set",title="LFR+RW - RW samples: $nsample")
#     histogram!(ls,bins=0:10:100,ylabel="",xlabel="size of planted set",label="",title="Total Sets: $(length(ls))",
#         xlims=(0,100),ylims=(0,500),
#         titlefontsize=10,
#         inset=bbox(0.7, 0.21, 0.2, 0.2, :bottom),
#         subplot=2)
#     Plots.savefig(p1,"purgatory/$(gname[1:end-5]).png")
# end


## size of symmetric difference in sets
# nsample = nsamples[1]
# gname = g*"-connect-graph-invdegdepth-"*string(nsample)*"-0.5-0.0-100-5.smat"
# A = loadGraph(gname,gpath)
# cmap = read_cmap(gpath*"cmap-gc-$(gname[1:end-5]).txt")
# mqimap = read_cmap(gpath*"cmap-gc-mqi-$(gname[1:end-5]).txt")
# slmap = read_cmap(gpath*"cmap-gc-sl-$(gname[1:end-5]).txt")
#
# sets = build_sets(cmap)
# mqisets = build_sets(mqimap)
# slsets = build_sets(slmap)
#
# ks = sort(collect(keys(sets)))
#
# rowids = rowvals(A)
# # c1 = map(x->conductance(A,sets[x],rowids=rowids),ks)
# # c2 = map(x->conductance(A,mqisets[x],rowids=rowids),ks)
# c3 = map(x->conductance(A,slsets[x],rowids=rowids),ks)
# # inds = findall(all.(zip(abs.(c1.-c3).>=0.3,c3.<0.1)))
# inds = findall(c3.<=0.05)
# ls = map(x->length(sets[x]),ks[inds])
# # delta = map(x->length(intersect(sets[x],slsets[x]))/length(slsets[x]),ks[inds])
# delta = map(x->length(intersect(sets[x],slsets[x]))/length(slsets[x]),ks[inds])
# delta = map(x->length(intersect(sets[x],slsets[x]))/length(sets[x]),ks)
# histogram(delta,leg=false,
#     xlabel="|(Planted) ⋂ (New)|/|Planted|",
#     ylabel="normalized frequency",
#     title="Planted vs SimpleLocal",
#     normalize=:probability,
#     ylims=(-0.01,0.5))
# Plots.savefig("purgatory/size-dist-sl-$g-1.png")
#
#
# p1 = scatter(inds,c1[inds],yscale=:log10,ylims=(1e-6,1.1),label="planted set",legend=:bottomleft,c=:black,xlims=(0,length(ks)+1))
# # scatter!(p1,inds,c2[inds],label="MQI")
# scatter!(p1,inds,c3[inds],label="SimpleLocal")
# plot!(p1,ylabel="conductance",xlabel="planted set",title="LFR+RW - RW samples: $nsample")
# histogram!(ls,bins=0:10:100,ylabel="",xlabel="size of planted set",label="",title="Total Sets: $(length(ls))",
#     xlims=(0,100),ylims=(0,500),
#     titlefontsize=10,
#     inset=bbox(0.7, 0.21, 0.2, 0.2, :bottom),
#     subplot=2)



#static images
# for numRWs = 500:500:8000
#     gname = g*"-connect-graph-invdegdepth-"*string(numRWs)*"-0.5-0.0-100-5.smat"
#     A = loadGraph(gname,gpath)
#     cmap = read_cmap(gpath*"cmap-gc-$(gname[1:end-5]).txt")
#     mqimap = read_cmap(gpath*"cmap-gc-mqi-$(gname[1:end-5]).txt")
#     slmap = read_cmap(gpath*"cmap-gc-sl-$(gname[1:end-5]).txt")
#
#     sets = build_sets(cmap)
#     mqisets = build_sets(mqimap)
#     slsets = build_sets(slmap)
#
#     ks = sort(collect(keys(sets)))
#
#     rowids = rowvals(A)
#     c1 = map(x->conductance(A,sets[x],rowids=rowids),ks)
#     c2 = map(x->conductance(A,mqisets[x],rowids=rowids),ks)
#     c3 = map(x->conductance(A,slsets[x],rowids=rowids),ks)
#     # inds = findall(all.(zip(abs.(c1.-c3).>=0.3,c3.<0.1)))
#     inds = findall(c3.<=0.05)
#     ls = map(x->length(sets[x]),ks[inds])
#     # delta = map(x->length(symdiff(sets[x],slsets[x])),ks[inds])
#     delta = map(x->length(intersect(sets[x],slsets[x]))/length(sets[x]),ks[inds])
#
#     p1 = scatter(inds,c1[inds],yscale=:log10,ylims=(1e-6,1.1),label="planted set",legend=:bottomleft,c=:black,xlims=(0,length(ks)+1))
#     # scatter!(p1,inds,c2[inds],label="MQI")
#     scatter!(p1,inds,c3[inds],label="SimpleLocal")
#     plot!(p1,ylabel="Conductance",
#         xlabel="Set",
#         title="LFR+RW - RW samples: $numRWs\nFinal Conductance ≤ 0.05\nTotal Sets: $(length(inds))")
#     histogram!(ls,bins=0:10:100,
#         normalize=:probability,
#         ylabel="",
#         xlabel="|Planted|",
#         label="",
#         xlims=(0,100),ylims=(-0.01,0.4),
#         xguidefontsize=10,
#         inset=bbox(0.4, 0.21, 0.2, 0.2, :bottom),
#         subplot=2)
#     histogram!(delta,bins=0:0.1:1.0,ylabel="",
#         xlabel="|Planted ⋂ New|/|Planted|",
#         label="",
#         normalize=:probability,
#         # title="Planted vs SimpleLocal Sets",
#         xlims=(0,1.0),ylims=(-0.01,0.4),
#         xguidefontsize=10,
#         inset=bbox(0.72, 0.215, 0.2, 0.2, :bottom),
#         subplot=3,
#         c=:green)
#
#     Plots.savefig("purgatory/planted-sets-$g-$numRWs.png")
# end
#
#
#
#
# using Interact
# wdg = slider(1:100)
# vbox(wdg, observe(wdg))


##sets missed in diffusion
# include("fast-diffusion.jl")
# include("ncpplots.jl")
# using DataFrames
#=
hexhist = fit(HexBinPlots.HexHistogram,
    log10.(szs),log10.(conds),100,density=true)
setinds = xy2inds_(log10.(szs),log10.(conds),100)

h,vh = HexBinPlots.make_shapes(hexhist)
xsh = Vector{Float64}()
ysh = Vector{Float64}()
for k in eachindex(h)
  append!(xsh,h[k].x)
  push!(xsh,h[k].x[1])
  push!(xsh,NaN)
  append!(ysh,h[k].y)
  push!(ysh,h[k].y[1])
  push!(ysh,NaN)
end

f = Plots.plot(10.0.^xsh[1:8],10.0.^ysh[1:8],fill_z=vh[1],linecolor=nothing,
    seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=true,
    ylims=(1e-2,1),xlims=(3,1.2*size(A,1)),leg=false)
for k=2:length(vh)
    inds = (k-1)*8+1:8*k
    Plots.plot!(10.0.^xsh[inds],10.0.^ysh[inds],fill_z=vh[k],linecolor=nothing,
        seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=true,
        ylims=(1e-2,1),xlims=(3,1.2*size(A,1)),leg=false)
end
f
ylims!(f,(1e-4,1))
=#

#we have a problem here. most of the sets we generate in
#the ncp have small size (5-50 nodes). how do i account for
#the bin imbalance when generating heatmap?

#=
should probably look at each bin.
for each bin, need some way to weigh sets in bin
in regards to missed nodes.

for each bin
    aggregate sets in that bin into Sbin
    weight that bin by how many nodes are susceptible in Sbin
    so get 0-1 measure for each bin and plot as heatmap

how to figure out what bin a set is in? adapt code from HexBinPlots
=#
#hexhistogram testing
# include("hexbins.jl")
#
# using Hexagons
# function xy2inds_(x::AbstractArray,y::AbstractArray,
#     bins::Union{NTuple{1,Int},NTuple{2,Int},Int})
#
#     if length(bins) == 2
#         xbins, ybins = bins
#     else
#         xbins, ybins = (bins...,bins...)
#     end
#     xmin, xmax = extrema(x)
#     ymin, ymax = extrema(y)
#     xspan, yspan = xmax - xmin, ymax - ymin
#     xsize, ysize = xspan / xbins, yspan / ybins
#     x0, y0 = xmin - xspan / 2,ymin - yspan / 2
#
#     inds = Dict{(Tuple{Int, Int}), Vector{Int}}()
#     cnts = Dict{(Tuple{Int, Int}), Int}()
#
#     for i in eachindex(x)
#         h = convert(HexagonOffsetOddR,
#                     cube_round(x[i] - x0 + xsize, y[i] - y0 + ysize, xsize, ysize))
#         idx = (h.q, h.r)
#
#         cnts[idx] = 1 + get(cnts,idx,0)
#         inds[idx] = push!(get(inds,idx,Vector{Int}()),i)
#     end
#     inds,xsize,ysize,x0,y0
# end
# function inds2xy_(inds::Dict{S,T},xsize, ysize, x0, y0) where {S,T}
#     nhex = length(inds)
#     xh = zeros(nhex)
#     yh = zeros(nhex)
#     vh = zeros(Int,nhex)
#     k = 0
#     for (idx, s) in inds
#         k += 1
#         xx,yy = Hexagons.center(HexagonOffsetOddR(idx[1], idx[2]),
#         xsize, ysize, x0, y0)
#         xh[k] = xx
#         yh[k] = yy
#         vh[k] = length(s)
#     end
#     xh,yh,vh
# end 
#
# function plot_missed_sets(gname::String,ncploc::String="",sloc::String="";
#             nbins::Int=100,ntrials::Int=1)
#     if endswith(gname,".smat")
#         g = gname[1:end-5]
#     else
#         g = gname
#     end
#     ncp,headerinfo = readdlm(ncploc*"ncpinfo-$g.txt",',',header=true)
#     ncp = DataFrame(ncp,vec(headerinfo))
#
#     ncpmat = MatrixNetworks.readSMAT(ncploc*"ncpmat-$g.smat")
#
#     #load snodes data
#     snodes = readdlm(sloc*"scounts-sir-0.05-0.05.txt",',',Int)
#     if size(snodes,2)>1
#         snodes = snodes[:,1] #using the first column #no quarantining
#     end
#     return plot_missed_sets(ncp,ncpmat,snodes;nbins=nbins,ntrials=ntrials)
# end
#
# function plot_missed_sets(ncp::DataFrame,ncpmat::SparseMatrixCSC,snodes;
#             nbins::Int=100,ntrials::Int=1)
#
#     @assert(maximum(snodes)<=ntrials,"ntrials < max # of times a node is susceptible")
#     #larger values means we missed more nodes in that set
#     missedsets = ncpmat*snodes./sum(ncpmat;dims=2)
#     #for aggregating over different diffusions
#     missedsets./=ntrials
#
#     x,y = log10.(ncp.size),log10.(ncp.cond)
#     inds,xsize,ysize,x0,y0 = xy2inds_(x,y,nbins)
#     xh,yh,vh = inds2xy_(inds,xsize,ysize,x0,y0)
#     hexhist = HexBinPlots.HexHistogram(xh,yh,vh,xsize,ysize,false)
#     h,vh = HexBinPlots.make_shapes(hexhist)
#
#     xsh = Vector{Float64}()
#     ysh = Vector{Float64}()
#     for k in eachindex(h)
#       append!(xsh,h[k].x)
#       push!(xsh,h[k].x[1])
#       push!(xsh,NaN)
#       append!(ysh,h[k].y)
#       push!(ysh,h[k].y[1])
#       push!(ysh,NaN)
#     end
#
#     #normalizing bins
#     wh = Vector{Float64}()
#     sets = ncpmat_to_sets(ncpmat)
#     for (idx,s) in inds
#         #fraction of susceptible nodes in bin
#         nodes = collect(union(sets[s]...))
#         push!(wh,sum(snodes[nodes])/length(nodes))
#
#         #mean_{sets in bin}(frac of set susceptible)
#         # push!(wh,sum(missedsets[s])/length(s))
#     end
#
#     #couldn't get fill_z to work. might have undergone a change since v1.5.1
#     #so will jerry-rig this to plot one by one.
#     f = Plots.plot(10.0.^xsh[1:8],10.0.^ysh[1:8],fill_z=wh[1],linecolor=nothing,
#         seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=true,
#         ylims=(1e-2,1),xlims=(3,1.2*size(A,1)),leg=false,
#         clims=(0,1))
#     for k=2:length(wh)
#         winds = (k-1)*8+1:8*k
#         Plots.plot!(10.0.^xsh[winds],10.0.^ysh[winds],fill_z=wh[k],linecolor=nothing,
#             seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=true,
#             ylims=(1e-2,1),xlims=(3,1.2*size(A,1)),leg=false,
#             clims=(0,1))
#     end
#     ylims!(f,(1e-4,1))
#     return f
# end
#
# function plot_missed_sets1(ncp::DataFrame,ncpmat::SparseMatrixCSC,snodes;
#             nbins::Int=100,ntrials::Int=1)
#
#     @assert(maximum(snodes)<=ntrials,"ntrials < max # of times a node is susceptible")
#     #larger values means we missed more nodes in that set
#     missedsets = ncpmat*snodes./sum(ncpmat;dims=2)
#     #for aggregating over different diffusions
#     missedsets./=ntrials
#
#     x,y = log10.(ncp.size),log10.(ncp.cond)
#     inds,xsize,ysize,x0,y0 = xy2inds_(x,y,nbins)
#     xh,yh,vh = inds2xy_(inds,xsize,ysize,x0,y0)
#     hexhist = HexBinPlots.HexHistogram(xh,yh,vh,xsize,ysize,false)
#     h,vh = HexBinPlots.make_shapes(hexhist)
#
#     xsh = Vector{Float64}()
#     ysh = Vector{Float64}()
#     for k in eachindex(h)
#       append!(xsh,h[k].x)
#       push!(xsh,h[k].x[1])
#       push!(xsh,NaN)
#       append!(ysh,h[k].y)
#       push!(ysh,h[k].y[1])
#       push!(ysh,NaN)
#     end
#
#     #normalizing bins
#     wh = Vector{Float64}()
#     # sets = ncpmat_to_sets(ncpmat)
#     for (idx,s) in inds
#         #fraction of susceptible nodes in bin
#         # nodes = collect(union(sets[s]...))
#         # push!(wh,sum(snodes[nodes])/length(nodes))
#
#         #mean_{sets in bin}(frac of set susceptible)
#         push!(wh,sum(missedsets[s])/length(s))
#     end
#
#     #couldn't get fill_z to work. might have undergone a change since v1.5.1
#     #so will jerry-rig this to plot one by one.
#     f = Plots.plot(10.0.^xsh[1:8],10.0.^ysh[1:8],fill_z=wh[1],linecolor=nothing,
#         seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=true,
#         ylims=(1e-2,1),xlims=(3,1.2*size(A,1)),leg=false,
#         clims=(0,1))
#     for k=2:length(wh)
#         winds = (k-1)*8+1:8*k
#         Plots.plot!(10.0.^xsh[winds],10.0.^ysh[winds],fill_z=wh[k],linecolor=nothing,
#             seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=true,
#             ylims=(1e-2,1),xlims=(3,1.2*size(A,1)),leg=false,
#             clims=(0,1))
#     end
#     ylims!(f,(1e-4,1))
#     return f
# end
#
#
# # gs = ["wiki-talk","brightkite","gowalla"]
# gs = ["mexico","dblp","anon"]
# gname = getgnames(gs[3],"graphs/real-graphs/")[1]
# ncp,headerinfo = readdlm("ncpdata/ncpinfo-$(gname[1:end-5]).txt",',',header=true)
# ncp = DataFrame(ncp,vec(headerinfo))
# ncpmat = MatrixNetworks.readSMAT("ncpdata/ncpmat-$gname")
#
# A = loadGraph("graphs/real-graphs/$gname")
# ncp.size = map(x->min(x,size(A,2)-x),ncp.size)

# beta in [0.01 0.05 0.1]

# beta = 0.01
# println("working on beta=$beta")
# E = EventEpidemic.SIRData(A;beta=beta)
# seed = rand(1:size(A,1))
# l,E = EventEpidemic.epidemic(E,seed)
# netinfs,newinfs = EventEpidemic.get_infdata(l,E)
# while sum(newinfs)<0.1*size(A,1)
#     l,E = EventEpidemic.epidemic(E,seed)
#     netinfs,newinfs = EventEpidemic.get_infdata(l,E)
# end
# sum(E.snodes)/length(E.snodes)
#
# f = plot_missed_sets(ncp,ncpmat,E.snodes)
# xlabel!(f,"size")
# ylabel!(f,"conductance")
# title!(f,"$(gname[1:end-5]) - SIR($beta,0.05)\nFraction of Infected Nodes: $(round(1-sum(E.snodes)/length(E.snodes),digits=2)*100)%")
# Plots.savefig(f,"purgatory/missed-sets/missed-sets-$(gname[1:end-5])-sir-$beta-0.05-node-metric.png")
#
# f = plot_missed_sets1(ncp,ncpmat,E.snodes)
# xlabel!(f,"size")
# ylabel!(f,"conductance")
# title!(f,"$(gname[1:end-5]) - SIR($beta,0.05)\nFraction of Infected Nodes: $(round(1-sum(E.snodes)/length(E.snodes),digits=2)*100)%")
# Plots.savefig(f,"purgatory/missed-sets/missed-sets-$(gname[1:end-5])-sir-$beta-0.05-set-metric.png")

# println("Performing aggreagted trials")
# beta = 0.1
# E = EventEpidemic.SIRData(A;beta=beta)
# scount = zeros(Int,length(E.snodes))
# k = 500
# for i = 1:k
#     println("working on diffusion $i")
#     EventEpidemic.rtrials(1,E,rand(1:size(A,1)),scount=scount,inflb=1)
# end
# sum(scount)/(k*length(E.snodes))
#
#
# f = plot_missed_sets(ncp,ncpmat,scount./k)
# xlabel!(f,"size")
# ylabel!(f,"conductance")
# title!(f,"$(gname[1:end-5]) - SIR($beta,0.05)\nAvg. Infected Fraction: $(round(1-sum(scount)/(k*length(E.snodes)),digits=2)*100)%")
# Plots.savefig(f,"purgatory/missed-sets/missed-sets-$(gname[1:end-5])-sir-$beta-0.05-node-metric-$k.png")

#aggregating over sets. aggregating over nodes looks very similar and is easier to interpret.
# f = plot_missed_sets1(ncp,ncpmat,scount./k)
# xlabel!(f,"size")
# ylabel!(f,"conductance")
# title!(f,"$(gname[1:end-5]) - SIR($beta,0.05)\nAvg. Infected Fraction: $(round(1-sum(scount)/(k*length(E.snodes)),digits=2)*100)%")
# Plots.savefig(f,"purgatory/missed-sets/missed-sets-$(gname[1:end-5])-sir-$beta-0.05-set-metric-$k.png")



##
# gpath = "graphs/real-graphs/"
# gnames = getgnames("soc",gpath)
# gname = gnames[5]
# A = loadGraph(gname,gpath)
# A = largest_component(A)[1]
#
# davg = nnz(A)/size(A,1)
# ncp,sets = make_ncp(A,get_sets=true)
#
#
# f = plot_hexbinncp(ncp,A)
# title!(f,"NCP - $(gname[1:end-5])\ndavg: $(round(davg,digits=2)), nnodes: $(size(A,1))")
# # Plots.savefig(f,dst*"ncphexbin-$(gname[1:end-5]).png")
#
#
#
# f = ncpcoverage(ncp,sets,A)
# title!(f,"Coverage - $(gname[1:end-5])\ndavg: $(round(davg,digits=2)), nnodes: $(size(A,1))")
# # Plots.savefig(f,dst*"ncpcoverage-$(gname[1:end-5]).png")
#
#
# #make ncps
# # ia-wiki-talk
# # brightkite
# # gowalla
# gs = ["ia-wiki-talk","brightkite","gowalla"]
# gname = getgnames(gs[3],gpath)[1]
# A = loadGraph(gname,gpath)
# A = largest_component(A)[1]
#
# davg = nnz(A)/size(A,1)
# ncp,sets = make_ncp(A,get_sets=true)
#
# #save raw data
# CSV.write("ncpdata/ncpinfo-$(gname[1:end-5]).txt",ncp)
# writeSMAT(makeSparseMatrix(sets,size(A,1)),"ncpdata/ncpmat-$gname")
# A
# ##
# include("fast-diffusion.jl")
# using Plots
# gpath = "graphs/real-graphs/"
# gname = getgnames("mexico",gpath)[1]
#
# function secondary_infs_node(l::Int,E,seed::Int)
#     """
#         uses infection history from E to count number of secondary infections for each node
#         ties are broken by counting all secondary infs for the node with the smallest label.
#         result[v] = nodes directly infected by v
#     """
#     indshist = first.(allmins.(E.ihist))
#     hist = map(x->E.inNeighborsA[x][indshist[x]],1:length(E.snodes))
#     # hist[v] = nodes that infect v if inf gets passed
#     result = Vector{Vector{Int}}()
#     for i=1:length(E.snodes)
#         push!(result,Vector{Int}())
#     end
#
#     @inbounds for v = 1:length(E.snodes)
#         if E.itime[v]<l #got infected
#             push!(result[hist[v][1]],v) #smallest labeled node that infs v
#         end
#     end
#     #fix seed node
#     deleteat!(result[hist[seed][1]],findfirst(result[hist[seed][1]].==seed))
#     return result
# end
#
# #use r0 and itime to compute rt
# function get_rt(lastactive,E,seed)
#     """
#         for time t=0:lastactive, compute the effective reproduction number
#         this is the quantity |∂(Iₜ)|/|Iₜ|
#     """
#
#     r0 = secondary_infs_node(lastactive,E,seed)
#     result = zeros(Float64,lastactive+1) #time t=0 gives r0 for seed node.
#
#     @inbounds for v = 1:length(E.snodes)
#         t1 = E.itime[v]
#         for u in r0[v]
#             t2 = E.itime[u]
#             result[t1+1:t2].+=1 #indices start at 0 so shift by 1
#         end
#     end
#
#     cuminfs = cumsum(EventEpidemic.get_infdata(lastactive,E)[2]).+1
#     result[1:end-1]./=cuminfs
#     return result
# end
#
# function vec_vec_to_array(input::Vector{Vector{T}}) where {T}
#     lens = length.(input)
#     maxlength = maximum(lens)
#     ms = zeros(Float64,maxlength)
#     stds = zeros(Float64,maxlength)
#     result = zeros(Float64,length(input),maxlength)
#
#     for t = 1:maxlength
#         for k = 1:length(input)
#             if lens[k]>=t
#                 result[k,t] = input[k][t]
#             end
#         end
#     end
#     return result
# end
#
# function make_rt_plots(gname,qpercent=0;beta::Float64=0.05,gpath::String="")
#     f = Plots.plot()
#     gs = [gname,"rewired-10000.0-$gname","er-10000.0-$gname"]
#     labdict = Dict([(1,"original"),(2,"rewired"),(3,"er")])
#     A = loadGraph("$gpath$gname")
#     qcap = round(Int,qpercent/100*size(A,1))
#
#     for (i,g) in enumerate(gs)
#         println("working on graph $g")
#         A = loadGraph("$gpath$g")
#         d = vec(sum(A;dims=2))
#
#         E = EventEpidemic.SIRData(A,beta=beta,qcap=qcap)
#         res = Vector{Vector{Float64}}()
#         for i = 1:500
#             seed = rand(1:size(A,1))
#             l,E = EventEpidemic.epidemic(E,seed,q=true)
#             push!(res,get_rt(l,E,seed))
#         end
#
#         rts = vec_vec_to_array(res)
#         ms = vec(mean(rts;dims=1))
#         stderrs = vec(mapslices(x->std(x)/sqrt(length(x)),rts;dims=1))
#         avginfs = mean(rts[:,end])/size(A,1)
#
#         # scatter!(f,1:length(ms),ms,yerr=stderrs,xscale = :log10,leg=false,
#         #         xlabel="time",ylabel="Effective Reproduction Number",
#         #         ylims=(-0.1,11),title="$(gname[1:end-5])\nSIR($beta,0.05)  nsamples: $(length(res))\nR0: ~$(round(ms[1],digits=2))")
#
#         scatter!(f,1:length(ms),ms,yerr=stderrs,xscale = :log10,leg=true,
#                 xlabel="time",ylabel="Effective Reproduction Number",
#                 ylims=(-0.1,11),title="$(gname[1:end-5])\nSIR($beta,0.05)  nsamples: $(length(res))\nAvg infs: $(avginfs)  Qcapacity: $qpercent%",
#                 label = labdict[i])
#     end
#     Plots.savefig(f,"purgatory/rt-plots/$(gname[1:end-5])-sir-$beta-0.05-$qpercent.png")
#     f
# end
#
# # gnames = ["modmexico","dblp","anon","enron"]
# # gnames = ["cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5",
# #     "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.9-0.0-100-5",
# #     "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17.smat"]
# # gnames = ["uf",
# #         "commutes","astro",
# gnames = ["penn","ucf",
#         "usf","modberkeley","modfsu","modharvard","modnotre-dame","modstanford",
#         "moduillinois","modunc","modwisconsin"]
# qpercents = [0,10]
# # gname = getgnames("mexico",gpath)[1]
# # f = make_rt_plots(gname,0,beta=0.05,gpath=gpath)
# ps = Vector()
# @showprogress for gname in gnames
#     println("working on $gname")
#     for qp in qpercents
#         gname = getgnames(gname,gpath)[1]
#         f = make_rt_plots(gname,qp,beta=0.05,gpath=gpath)
#         push!(ps,f)
#         if "cn-$gname" in getgnames(gname,gpath)
#             f = make_rt_plots("cn-$gname",qp,beta=0.05,gpath=gpath)
#             push!(ps,f)
#         end
#     end
# end
#
# os = deepcopy(ps)
# #adjust ylims
# for f in os
#     plot!(f,ylims=(-0.1,13))
# end
#
# #add in annotations
# for f in os
#     plot!(f,[1;1],collect(ylims(f)), annotations = ([1.3], [10.0], "R0"),
#         label="",linestyle=:dashdot,c=:black)
#     plot!(f,collect(xlims(f)),[1;1], annotations = ([300], [1.5], "Rt=1"),
#         label="",linestyle=:dashdot,c=:black)
# end
# os[1]
# os[2]
# os[3]
# os[4]
# os[50]
#
# for (i,val) in enumerate(vcat(Base.Iterators.product(qpercents,gnames)...))
#     Plots.savefig(os[i],"purgatory/rt-plots/$(val[2])-sir-0.05-0.05-$(val[1]).png")
# end




##making table for graph statistics/features/attributes
# using MatrixNetworks
# using LightGraphs
# gpath = "graphs/real-graphs/"
# gs = ["enron","uf","modmexico","anony",
#         "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5",
#         "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.9-0.0-100-5",
#         "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17.smat",
#         "dblp","commutes","astro","penn","ucf",
#         "usf","modberkeley","modfsu","modharvard","modnotre-dame","modstanford",
#         "moduillinois","modunc","modwisconsin","flickr","livejournal"]

#only keep original graph and sparsified graph. maybe do two rewired graphs or just the original graph

#just the Original + sparsified graphs
#gname,davg,
# res = Array{String}(["gname" "nnodes" "nedges" "davg" "d₂ avg" "average CC"])
# for g in gs
#     gnames = getgnames(g,gpath)
#     if startswith(g,"cl-")
#         filter!(x->!startswith(x,"er-") && !startswith(x,"rewired-") && !startswith(x,"kcore"),gnames)
#     else
#         filter!(x->!startswith(x,"er-") && !startswith(x,"rewired-") && !startswith(x,"cl-") && !startswith(x,"kcore"),gnames)
#     end
#     for gname in gnames
#         println("working on graph $(gname)")
#         A = loadGraph(gname,gpath)
#         d = Int.(sum(A;dims=2))
#
#         #found bug in clustercoeffs from MatrixNetworks for nodes of degree 1. sometimes gives non-zero cc for those nodes. see modStanford.
#         #use LightGraphs for that one instead
#         avgcc = mean(LightGraphs.local_clustering_coefficient(LightGraphs.Graph(A)))
#         res = vcat(res,[gname[1:end-5] "$(size(A,1))" "$(Int(nnz(A))/2)" "$(mean(d))" "$(mean(A*d))" "$(avgcc)"])
#     end
# end


#write to file
# writedlm("paper/graph-stats.csv",res,',')



function commonNeighbors(A::SparseMatrixCSC)
  @assert issymmetric(A)
  P = deepcopy(A)
  fill!(nonzeros(P),0)
  rowval = rowvals(A)
  neighs = zeros(Bool, size(A,1))
  for j=1:size(A,1)
    # index neighbors
    for nzj in nzrange(A, j)
      neighs[rowval[nzj]] = true
    end
    # for each edge out of this node...
    for nzj in nzrange(A, j)
      i = rowval[nzj]
      score = 0
      for nzi in nzrange(A, i)
        w = rowval[nzi]
        if neighs[w] == true
          # for each neighbor of i (w) that is
          # also a neighbor of j neighs[w] = true
          # increase the score
          score += 1
        end
      end
      nonzeros(P)[nzj] = score
    end
    # reset indicies
    for nzj in nzrange(A, j)
      neighs[rowval[nzj]] = false
    end
  end
  return P
end

using Plots
using Measures

gpath = "homes/oeldagha/network-epidemics/data/graphs/"
dloc="data/test/high-percentage-transfer/seir/"
ntrials = 50
# qpercents = [0;5;10;15;20;25;collect(30:10:100)]
method="seir"
function graph_qpercent_heatmap(gname::String,betas::Vector{Float64},gammas::Vector{Float64}=[0.05];
                        gpath::String="",dtype::String="tinfs",dloc::String="",
                        method::String="sir",
                        rps::Vector{Float64}=vec([0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.25 0.5 1.0 10.0]),
                        ntrials::Int=500,
                        qpercents::Vector{Int}=collect(0:15))

    gr()
    # do stuff
    A = loadGraph(gname,gpath)
    nnodes = size(A,1)
    ths = nnodes.*collect(0.01:0.01:1.0)
    g = gname[1:end-5]
    qfigs,figs = Vector(),Vector()
    for beta in betas
        for gamma in gammas
            data = Vector{Vector{Vector{Int}}}()
            gnames = Vector{String}()
            for p in 100 .*reverse(rps)# 0.1 0.25 0.5 1.0]# 10.0]#[1 5 10 25 50 100 1000]
                push!(data,cumsum.(read_inf_data("rewired-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype,method=method)))
                push!(gnames,"rewired-$p-$g")
            end
            push!(data,cumsum.(read_inf_data(g,dloc=dloc,beta=beta,gamma=gamma,dtype=dtype,method=method)))
            push!(gnames,g)
            for p in 100 .*rps# 0.1 0.25 0.5 1.0]# 10.0]#[1 5 10 25 50 100 1000]
                push!(data,cumsum.(read_inf_data("er-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype,method=method)))
                push!(gnames,"er-$p-rewired-$g")
            end

            #build heatmap
            #will do color as pr( >=k% of nodes)
            pdata = zeros(length(ths),length(gnames))
            for graphind=1:length(gnames)
                if dtype == "tinfs"
                    ps = map(x->findlast(x.>ths),last.(data[graphind][1:ntrials]))
                elseif dtype == "cinfs"
                    ps = map(x->findlast(x.>ths),maximum.(data[graphind][1:ntrials]))
                end
                ps = Vector{Union{Int,Nothing}}(ps)
                ps[ps.==nothing].= 0
                ps = Int.(ps)
                pdata[:,graphind] .= fit(Histogram,ps,1:length(ths)+1).weights
                pdata[:,graphind] .= cumsum(pdata[end:-1:1,graphind])[end:-1:1]./ntrials #500 was the sample size for the data
            end
            if dtype == "tinfs"
                ylab = "Pr(Total Infs ≥ k% of Nodes)"
            elseif dtype == "cinfs"
                ylab = "Pr(Maximum Infs ≥ k% of Nodes)"
            end
            f1 = heatmap(pdata,ylabel="Percent of Total Nodes",colorbar_title=ylab,
                        title="$(gname[1:end-5]) - $(uppercase(method))($beta,$gamma)",clims=(0,1),
                        xlabel="     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent",
                        xrotation = 90)

            if dtype == "tinfs"
                clabel = "Fraction of Infected Nodes"
                data = hcat(map(x->last.(x),data)...)./nnodes
            elseif dtype == "cinfs"
                clabel = "Normalized Maximum Infection Size"
                data = hcat(map(x->maximum.(x),data)...)./nnodes
            end
            f = heatmap(data,ylabel="quarantine percentage",
                xlabel="     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent",
                title="$(gname[1:end-5]) - $(uppercase(method))($beta,$gamma)\n$(ntrials) node sample",
                colorbar_title=clabel,c=:heat,clims=(0,1),
                yticks=(collect(round(Int,ntrials/2):ntrials:size(data,1)),qpercents),xrotation = 90)

            plot!(f,xticks=(1:2*length(rps)+1, string.(round.(100 .*vec([reverse(rps)... 0 rps...]),digits=2))),margins=9Measures.mm)
            plot!(f1,xticks=(1:2*length(rps)+1, string.(round.(100 .*vec([reverse(rps)... 0 rps...]),digits=2))),margins=9Measures.mm)
            push!(qfigs,f)
            push!(figs,f1)
        end
    end
    # plot!(f,xticks=(1:2*length(rps), string.(round.(100 .*vec([reverse(rps)... rps...]),digits=2))),margins=9mm)
    # plot!(f1,xticks=(1:2*length(rps), string.(round.(100 .*vec([reverse(rps)... rps...]),digits=2))),margins=9mm)
    return qfigs,figs
end

function graph_qpercent_heatmap1(gname::String,betas::Vector{Float64},gammas::Vector{Float64}=[0.05];
    gpath::String="",dtype::String="tinfs",dloc::String="",
    method::String="sir",
    rps::Vector{Float64}=vec([0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.25 0.5 1.0 10.0]),
    ntrials::Int=50,
    qpercents::Vector{Int}=collect(0:15))

    # gr()
    # do stuff
    A = loadGraph(gname,gpath)
    nnodes = size(A,1)
    ths = nnodes.*collect(0.01:0.01:1.0)
    g = gname[1:end-5]
    qfigs = Vector()
    for beta in betas
        for gamma in gammas
            data = Vector{Vector{Vector{Int}}}()
            gnames = Vector{String}()
            for p in 100 .*reverse(rps)# 0.1 0.25 0.5 1.0]# 10.0]#[1 5 10 25 50 100 1000]
                push!(data,cumsum.(read_inf_data("rewired-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype,method=method)))
                push!(gnames,"rewired-$p-$g")
            end
            push!(data,cumsum.(read_inf_data(g,dloc=dloc,beta=beta,gamma=gamma,dtype=dtype,method=method)))
            push!(gnames,g)
            for p in 100 .*rps# 0.1 0.25 0.5 1.0]# 10.0]#[1 5 10 25 50 100 1000]
                push!(data,cumsum.(read_inf_data("er-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype,method=method)))
                push!(gnames,"er-$p-rewired-$g")
            end
            if dtype == "tinfs"
                clabel = "Fraction of Infected Nodes"
                data = hcat(map(x->last.(x),data)...)./nnodes
            elseif dtype == "cinfs"
                clabel = "Normalized Maximum Infection Size"
                data = hcat(map(x->maximum.(x),data)...)./nnodes
            end
            #aggregating over nodes in a cell. 
            # println(size(data))
            # println(size(data))
            tmp = zeros(Float64,Int(size(data,1)/ntrials),size(data,2))
            @inbounds for col = 1:size(tmp,1)
                for row =1:size(tmp,1)
                    offset = (row-1)*ntrials
                    tmp[row,col] = sum(data[offset+1:offset+ntrials,col]/ntrials)
                end
            end
            f = heatmap(tmp,ylabel="quarantine percentage",
            xlabel="     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent",
            title="$(gname[1:end-5]) - $(uppercase(method))($beta,$gamma)\n$(ntrials) node sample",
            colorbar_title=clabel,c=:heat,clims=(0,1),
            yticks=(0.5:length(qpercents)-0.5,qpercents),xrotation = 90)

            plot!(f,xticks=(1:2*length(rps)+1, string.(round.(100 .*vec([reverse(rps)... 0 rps...]),digits=2))),margins=9Measures.mm)
            push!(qfigs,f)
        end
    end
    return qfigs
end


function graph_qpercent_heatmap2(gname::String,beta::Float64=0.05,gamma::Float64=0.05;
    gpath::String="",dtype::String="tinfs",dloc::String="",
    method::String="sir",
    rps::Vector{Float64}=vec([0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.25 0.5 1.0 10.0]),
    ntrials::Int=50,
    qpercents::Vector{Int}=collect(0:15))

    pyplot()
    # do stuff
    A = loadGraph(gname,gpath)
    nnodes = size(A,1)
    ths = nnodes.*collect(0.01:0.01:1.0)
    g = gname[1:end-5]
   
    data = Vector{Vector{Vector{Int}}}()
    gnames = Vector{String}()
    for p in 100 .*reverse(rps)# 0.1 0.25 0.5 1.0]# 10.0]#[1 5 10 25 50 100 1000]
        push!(data,cumsum.(read_inf_data("rewired-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype,method=method)))
        push!(gnames,"rewired-$p-$g")
    end
    push!(data,cumsum.(read_inf_data(g,dloc=dloc,beta=beta,gamma=gamma,dtype=dtype,method=method)))
    push!(gnames,g)
    for p in 100 .*rps# 0.1 0.25 0.5 1.0]# 10.0]#[1 5 10 25 50 100 1000]
        push!(data,cumsum.(read_inf_data("er-$p-$g",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype,method=method)))
        push!(gnames,"er-$p-rewired-$g")
    end
    if dtype == "tinfs"
        clabel = "Fraction of Infected Nodes"
        data = hcat(map(x->last.(x),data)...)./nnodes
    elseif dtype == "cinfs"
        clabel = "Normalized Maximum Infection Size"
        data = hcat(map(x->maximum.(x),data)...)./nnodes
    end
    #aggregating over nodes in a cell. 
    # println(size(data))
    # println(size(data))
    #TODO should probably put this part in a separate fcn.
    nrows,ncols = size(data)
    nrows = Int(nrows/ntrials)
    tmp = zeros(Float64,nrows,ncols)    
    for col = 1:ncols
        for row = 1:nrows
            offset = (row-1)*ntrials
            tmp[row,col] = sum(data[offset+1:offset+ntrials,col])/ntrials
        end
    end
    f = heatmap(tmp,ylabel="quarantine percentage",
    xlabel="     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent",
    title="$(gname[1:end-5]) - $(uppercase(method))($beta,$gamma)\n$(ntrials) node sample",
    colorbar_title=clabel,c=:heat,clims=(0,1),
    yticks=(0.5:length(qpercents)-0.5,qpercents),xrotation = 90)

    plot!(f,xticks=(1:2*length(rps)+1, string.(round.(100 .*vec([reverse(rps)... 0 rps...]),digits=2))),margins=9Measures.mm)
    return data,tmp,f
end


#total steps 1e6
function graph_sparsification_heatmap(gname::String,betas::Vector{Float64},gammas::Vector{Float64}=[0.05];
                        gpath::String="",dtype::String="tinfs",dloc::String="",
                        method::String="sir",
                        rps::Vector{Float64}=collect(0.1:0.1:0.9),
                        ntrials::Int=500,
                        qpercents::Vector{Int}=collect(0:15))
    pyplot()
    # do stuff
    A = loadGraph(gname,gpath)
    nnodes = size(A,1)
    ths = nnodes.*collect(0.01:0.01:1.0)
    g = gname[1:end-5]

    qfigs,figs = Vector(),Vector()
    for beta in betas
        for gamma in gammas
            data = Vector{Vector{Vector{Int}}}()
            gnames = Vector{String}()
            push!(data,cumsum.(read_inf_data(g,dloc=dloc,beta=beta,gamma=gamma,dtype=dtype,method=method)))
            push!(gnames,g)
            for p in rps# 0.1 0.25 0.5 1.0]# 10.0]#[1 5 10 25 50 100 1000]
                push!(data,cumsum.(read_inf_data("$g-$p-7",dloc=dloc,beta=beta,gamma=gamma,dtype=dtype,method=method)))
                push!(gnames,"$g-$p-7")
            end

            #build heatmap
            #will do color as pr( >=k% of nodes)
            pdata = zeros(length(ths),length(gnames))
            for graphind=1:length(gnames)
                if dtype == "tinfs"
                    ps = map(x->findlast(x.>ths),last.(data[graphind][1:ntrials]))
                elseif dtype == "cinfs"
                    ps = map(x->findlast(x.>ths),maximum.(data[graphind][1:ntrials]))
                end
                ps = Vector{Union{Int,Nothing}}(ps)
                ps[ps.==nothing].= 0
                ps = Int.(ps)
                pdata[:,graphind] .= fit(Histogram,ps,1:length(ths)+1).weights
                pdata[:,graphind] .= cumsum(pdata[end:-1:1,graphind])[end:-1:1]./ntrials #500 was the sample size for the data
            end
            if dtype == "tinfs"
                ylab = "Pr(Total Infs ≥ k% of Nodes)"
            elseif dtype == "cinfs"
                ylab = "Pr(Maximum Infs ≥ k% of Nodes)"
            end
            f1 = heatmap(pdata,ylabel="Percent of Total Nodes",colorbar_title=ylab,
                        title="$(gname[1:end-5]) - $(uppercase(method))($beta,$gamma)",clims=(0,1),
                        xlabel="     Original ⇒ Rewired Graph\nRewiring Fraction - Total Steps: 1e6",
                        xrotation = 90)

            if dtype == "tinfs"
                clabel = "Fraction of Infected Nodes"
                data = hcat(map(x->last.(x),data)...)./nnodes
            elseif dtype == "cinfs"
                clabel = "Normalized Maximum Infection Size"
                data = hcat(map(x->maximum.(x),data)...)./nnodes
            end


            f = heatmap(data,ylabel="quarantine percentage",
                xlabel="     Original ⇒ Rewired Graph\nRewiring Fraction - Total Steps: 10^6",
                title="$(gname[1:end-5]) - $(uppercase(method))($beta,$gamma)\n$(ntrials) node sample",
                colorbar_title=clabel,c=:heat,clims=(0,1),
                yticks=(collect(round(Int,ntrials/2):ntrials:size(data,1)),qpercents),xrotation = 90)

            plot!(f,xticks=(1:length(rps)+1, string.([0;rps])),margins=9Measures.mm)
            plot!(f1,xticks=(1:length(rps)+1, string.([0;rps])),margins=9Measures.mm)

            push!(qfigs,f)
            push!(figs,f1)
        end
    end

    #build heatmap
    #will do color as pr( >=k% of nodes)
    qfigs,figs
end

include("fast-diffusion.jl")
gpath = "data/graphs/"
method="sir"
dloc = "data/infdata/test/sir/"
betas = [collect(1e-3:1e-3:1e-2);collect(2e-2:1e-2:1e-1)];
gname = getgnames("dblp",gpath)[1]

# pyplot()
d,td,f = graph_qpercent_heatmap2(gname,gpath=gpath,dloc=dloc,
                        method=method,ntrials=50)
f

f2 = contourf(1:23,1:16,(x,y)->td[y,x],clims=(0,1))
plot!(f2,colorbar_ticks=:auto,colorbar_title="test text")


using Plots
pyplot()
colorbar_ticks=(0:0.2:1.0, string.(round.(Int, (0:0.2:1.0) .* 100), "%"))
Plots.scatter(rand(10), rand(10), marker_z = rand(10), colorbar_ticks=colorbar_ticks,clim=(0,1), colorbar_title="My title")


