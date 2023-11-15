parent_dir = "/p/mnt/scratch/network-epi/"
include(joinpath(parent_dir,"code/graph-io.jl"))
include(joinpath(parent_dir,"code/fast-diffusion.jl"))
include(joinpath(parent_dir,"code/rt.jl"))

using SparseArrays
using DelimitedFiles
using Plots
using Measures
ENV["GKSwstype"] = "100"
#graph rendering function from david
function graphlines(A::SparseMatrixCSC,xy::AbstractArray{T,2}) where T
    px,py = zeros(T,0),zeros(T,0)
    P = [px,py]
    rows = rowvals(A)
    skip = NaN.*xy[:,begin] # first row
    for j=1:size(A,2) # for each column
        for nzi in nzrange(A, j)
            i = rows[nzi]
            if i > j
                push!.(P, @view xy[:,i])
                push!.(P, @view xy[:,j])
                push!.(P, skip)
            end
        end
    end
    return px, py
end

function epidemic_snapshot(A,xy,infection_times,t)
  f = Plots.plot(graphlines(A, xy)..., frame_style=:none,
    legend=false, linewidth=1, size=(200,200), linecolor=RGBA(0.5,0.5,0.5,0.15))
  scatter!(f,xy[1,:],xy[2,:], color=RGBA(0.5,0.5,0.5,1.0), 
    markerstrokewidth=0, markersize=1,
    frame_style=:none,margin=-20mm,
    legend=false, linewidth=1, size=(200,200),)
  #find relevant infected nodes
  infected = findall(infection_times.<t)
  #scatter 
  Plots.scatter!(f,xy[1,infected],xy[2,infected], 
    markerstrokewidth=0, markersize=2.5, color=2)
  plot!(f,dpi=300)
  return f
end

# Rt := average number of new infs caused by an infected individual at time t 
# so at time t, look at Iₜ and compute cut(Iₜ)/size(Iₜ) <- very similar to expansion. 
# see definition on top of page 13 from that siam review paper: https://arxiv.org/pdf/2004.09608.pdf


#load graph 
# gpath = joinpath(parent_dir,"input/graphs/")
# A = loadGraph("study-20-draft-150.smat",gpath)
# xy = readdlm(gpath*"study-20-draft-150.xy")

#=
#an intiial rendering of the graph 
# f = Plots.plot(graphlines(A, xy)..., frame_style=:none,
#   legend=false, linewidth=1, size=(200,200), linecolor=RGBA(0.4,0.4,0.4,0.15))
# Plots.scatter!(f,xy[1,:],xy[2,:], color=RGBA(0.1,0.1,0.1,1.0),margin=-20mm, 
#   markerstrokewidth=0, markersize=1)
pts = collect(zip(xy[1,:],xy[2,:]))
xpt,ypt = 0,1/2

ind = argmin(map(pt->(pt[1]-xpt)^2+(pt[2]-ypt)^2,pts))
# Plots.scatter!(f,xy[1,ind:ind],xy[2,ind:ind], color=RGBA(0.1,0.1,0.1,1.0), 
#   markerstrokewidth=0, markersize=3)

target = findall(map(pt->pt[1]<0.2 && 0.5<pt[2]<0.75,pts))
other = setdiff(1:size(A,1),target)

E = EventEpidemic.SEIRData(A,beta=0.015)
# @showprogress for rseed = 1:1000
#   Random.seed!(rseed)
#   l,E = EventEpidemic.epidemic(E,ind)
#   while sum((EventEpidemic.get_infdata(l,E))[2])<50
#     l,E = EventEpidemic.epidemic(E,ind)
#   end
#   #check that we first infect the target nodes 
#   if minimum(E.itime[other])>maximum(E.itime[target])
#     println("found random seed: rseed = $rseed")
#   end
# end
# minimum(E.itime[other])
# maximum(E.itime[target])

Random.seed!(915)
l,E = EventEpidemic.epidemic(E,ind)
while sum((EventEpidemic.get_infdata(l,E))[2])<50
    l,E = EventEpidemic.epidemic(E,ind)
end

#find nodes not in target set that get directly infected by a target node
r0 = secondary_infs_node(l,E,ind)
target_node_infs = unique(vcat(r0[target]...))
target_node_infs = setdiff(target_node_infs,target)
rt = get_rt(l,E.itime,ind)

#sanity check 
length(target_node_infs)/length(target)
rt[101]

f = Plots.plot(graphlines(A, xy)..., frame_style=:none,
    legend=false, linewidth=0.5, size=(200,200), linecolor=RGBA(0.5,0.5,0.5,0.15))
#add all nodes 
scatter!(f,xy[1,:],xy[2,:], color=RGBA(0.5,0.5,0.5,1.0), 
    markerstrokewidth=0, markersize=1,
    frame_style=:none,margin=-20mm,
    legend=false, linewidth=1, size=(200,200),dpi=300)
#add infected 
scatter!(f,xy[1,target],xy[2,target],
        markersize=2,markerstrokewidth=0,color=2)
#add all infected 
infected = findall(E.itime.<l)
scatter!(f,xy[1,infected],xy[2,infected],
        markersize=2,markerstrokewidth=0,color=2)
#add secondary infs from target set 
scatter!(f,xy[1,target_node_infs],xy[2,target_node_infs],
    markersize=4,markerstrokewidth=0,color=3)
#add secondary infs from target set as infected  
scatter!(f,xy[1,target_node_infs],xy[2,target_node_infs],
    markersize=4,markerstrokewidth=0,color=2)



### single plot 
rt = get_rt(l,E.itime,ind)
last_positive_ind = findlast(rt.>0)
plot(1:last_positive_ind,rt[1:last_positive_ind],leg=false,
    xlabel="time step (t)",ylabel="Rt",linestyle=:dash)
scatter!(1:last_positive_ind,rt[1:last_positive_ind],alpha=0.25)
scatter!(101:101,rt[101:101],c=:green,markersize=8,xscale=:log10,xlims=(0.9,1e3))

=#
### many samples
# include(joinpath(parent_dir,"code/ncp/ncpplots1.jl"))
# include(joinpath(parent_dir,"code/rewiring-functions.jl"))

# gname = getgnames("mexico","input/graphs/")[1]
# A = loadGraph(gname,"input/graphs/")

# #using different seed nodes
# using ProgressMeter
# function aggregated_rt(A::SparseMatrixCSC)
#     #set things up 
#     E = EventEpidemic.SEIRData(A,beta=0.08)

#     tstamps = Vector{Int}()
#     rts = Vector{Float64}()
#     @showprogress for k = 1:1000
#         seednode = rand(1:lastindex(E.snodes))
#         l,E = EventEpidemic.epidemic(E,seednode)
#         ktrials = 1
#         while sum((EventEpidemic.get_infdata(l,E))[2])<50 && ktrials<=10
#             l,E = EventEpidemic.epidemic(E,seednode)
#             ktrials += 1
#         end

#         if sum((EventEpidemic.get_infdata(l,E))[2])>=50
#             sample_rt = get_rt(l,E,seednode)
#             lastind = findlast(sample_rt.>0)
#             push!(rts,sample_rt[1:lastind]...)
#             push!(tstamps,collect(1:lastind)...)
#         end
#     end
#     return tstamps,rts 
# end

#TODO seed experiments and save rt data 

# ts,rts = aggregated_rt(A)
# f = myncpplot1(ts,rts,xlims=(0.9,1e3),ylims=(1e-6,2*maximum(rts)),plotmin=false)

# B = rewire_graph(A,50*nnz(A))
# ts,rts = aggregated_rt(B)
# f = myncpplot1(ts,rts,xlims=(0.9,1e3),ylims=(1e-3,2*maximum(rts)),plotmin=false)

# B = er_rewire_graph(A,50*nnz(A))
# ts,rts = aggregated_rt(B)
# f = myncpplot1(ts,rts,xlims=(0.9,1e3),ylims=(1e-3,2*maximum(rts)),plotmin=false)


# get_rt(l,E.itime,seednode)
# rt1 = get_rt(l,E.itime,seednode)
# lastind1 = findlast(rt1.>0)


# tstamps = Vector{Int}()
# rt = Vector{Float64}()
# @showprogress for k = 1:1000
#     seednode = rand(1:lastindex(E.snodes))
#     l,E = EventEpidemic.epidemic(E,seednode)
#     ktrials = 1
#     while sum((EventEpidemic.get_infdata(l,E))[2])<50 && ktrials<=10
#         l,E = EventEpidemic.epidemic(E,seednode)
#         ktrials += 1
#     end

#     if sum((EventEpidemic.get_infdata(l,E))[2])>50
#         sample_rt = get_rt(l,E.itime,seednode)
#         lastind = findlast(sample_rt.>0)
#         push!(rt,sample_rt[1:lastind]...)
#         push!(tstamps,collect(1:lastind)...)
#     end
# end

# f = myncpplot1(ts,rts,xlims=(0.9,1e3),ylims=(1e-3,2*maximum(rts)),plotmin=false)
# plot!(f,xlabel="time step (t)",ylabel="Rt")








# tstamps = Vector{Int}()
# rts = Vector{Float64}()
# @showprogress for k = 1:1000
#     seednode = rand(1:lastindex(E.snodes))
#     l,E = EventEpidemic.epidemic(E,seednode)
#     ktrials = 1
#     while sum((EventEpidemic.get_infdata(l,E))[2])<50 && ktrials<=10
#         l,E = EventEpidemic.epidemic(E,seednode)
#         ktrials += 1
#     end

#     if sum((EventEpidemic.get_infdata(l,E))[2])>50
#         sample_rt = get_rt(l,E.itime,seednode)
#         lastind = findlast(sample_rt.>0)
#         push!(rts,sample_rt[1:lastind]...)
#         push!(tstamps,collect(1:lastind)...)
#     end
# end




### Rt data 
using Distributed 
# addprocs(50)
parent_dir = "/p/mnt/scratch/network-epi/"
@everywhere @eval include(joinpath($parent_dir,"code/rt.jl"))

nprocs()
gnames = Vector{String}()
rt_param_dict = Dict{String,Float64}()
#gname -> (beta, qpercent) #got this from images 



for gname in getgnames("study","input/graphs/")
    rt_param_dict[gname]=0.03
end
for gname in getgnames("cl-lfr","input/graphs/")
    rt_param_dict[gname]=0.175
end
for gname in getgnames("cn-","input/graphs/")
    rt_param_dict[gname]=0.1
end
rt_param_dict[getgnames("anony","input/graphs/")[1]]=0.1
rt_param_dict[getgnames("dblp","input/graphs/")[1]]=0.08
rt_param_dict[getgnames("gowalla","input/graphs/")[1]]=0.08
rt_param_dict[getgnames("brightkite","input/graphs/")[1]]=0.05
rt_param_dict[getgnames("epinions","input/graphs/")[1]]=0.05
rt_param_dict[getgnames("flickr","input/graphs/")[1]]=0.05
rt_param_dict[getgnames("commutes","input/graphs/")[1]]=0.002
rt_param_dict[getgnames("mexico","input/graphs/")[1]]=0.05
rt_param_dict[getgnames("covidflows-2020_08_31.smat","input/graphs/")[1]]=0.0009
rt_param_dict[getgnames("covidflows-2020_08_31-filtered-20.smat","input/graphs/")[1]]=0.03
rt_param_dict[getgnames("geometric-100000-2d-1-20","input/graphs/")[1]]=0.05



@everywhere rt_param_dict = $rt_param_dict
#for each graph, add extremal rewired graphs 
gnames = Vector{String}()
for g in keys(rt_param_dict)
    gname = getgnames(g,"input/graphs/")[1]
    push!(gnames,gname)
    push!(gnames,"er-10000.0-$gname")
    push!(gnames,"rewired-10000.0-$gname")
end


## testing 
# gname = getgnames("study-11-2022-45","input/graphs/")[1]
# res = aggregated_rt(gname,beta=rt_param_dict[gname])

#gather data
res = @showprogress pmap(gname->aggregated_rt(gname,
                                    beta=rt_param_dict[canonical_graph_name(gname)],
                                    qpercent=15), 
                                    gnames)


 

# using Plots
# using Measures
# ENV["GKSwstype"] = "100"

# load data 
function load_rtdata(gname::String,beta::Float64=0.05)
    g = canonical_graph_name(gname)
    
    fpath = joinpath(parent_dir,"pipeline/data/$(g[1:end-5])/rt-1000-seir-$beta-0.05-$(gname[1:end-5]).txt")
    Rt = readdlm(fpath)

    fpath = joinpath(parent_dir,"pipeline/data/$(g[1:end-5])/rt-contagious-1000-seir-$beta-0.05-$(gname[1:end-5]).txt")
    contagious_Rt = readdlm(fpath)

    fpath = joinpath(parent_dir,"pipeline/data/$(g[1:end-5])/rt-exposed-1000-seir-$beta-0.05-$(gname[1:end-5]).txt")
    exposed_Rt = readdlm(fpath)

    return Int.(Rt[:,1]),Rt[:,2],contagious_Rt[:,2],exposed_Rt[:,2] #timestep, Rt, contagious_Rt, exposed_Rt
end

function load_rtdata(gname::String,beta::Float64=0.05,qpercent::Int)
    g = canonical_graph_name(gname)
    
    fpath = joinpath(parent_dir,"pipeline/data/$(g[1:end-5])/rt-1000-seir-$beta-0.05-$qpercent-$(gname[1:end-5]).txt")
    Rt = readdlm(fpath)

    fpath = joinpath(parent_dir,"pipeline/data/$(g[1:end-5])/rt-contagious-1000-seir-$beta-0.05-$qpercent-$(gname[1:end-5]).txt")
    contagious_Rt = readdlm(fpath)

    fpath = joinpath(parent_dir,"pipeline/data/$(g[1:end-5])/rt-exposed-1000-seir-$beta-0.05-$qpercent-$(gname[1:end-5]).txt")
    exposed_Rt = readdlm(fpath)

    return Int.(Rt[:,1]),Rt[:,2],contagious_Rt[:,2],exposed_Rt[:,2] #timestep, Rt, contagious_Rt, exposed_Rt
end


include("../ncp/ncpplots1.jl")

# gname = getgnames("mexico","input/graphs/")[1]

function rt_plots(gname::String)

    result = []
    
    if startswith(gname,"er-10000.0-")
        g = gname[12:end]
    elseif startswith(gname,"rewired-10000.0-")
        g = gname[17:end]
    else
        g = gname
    end
    
    fpath = joinpath(parent_dir,"pipeline/data/$(g[1:end-5])/rt-1000-seir-0.05-0.05-$(gname[1:end-5]).txt")
    data = readdlm(fpath)
    ts,rts = data[:,1],data[:,2]
    inds = findall(ts.==1)
    
    xmin,xmax = 0.9,305
    ymin,ymax = 1e-3,1e3

    # plot type 1 #vanilla plot w/ alpha thresholding
    f = Plots.plot(leg=false)
    for k = 1:min(100,lastindex(inds)-1)
        tmprt = rts[inds[k]:inds[k+1]-1]
        Plots.plot!(f,ts[inds[k]:inds[k+1]-1],tmprt,color=2,alpha=min(10/k,1),
            linestyle=:dash)
    end
    Plots.plot!(f,
        xscale=:log10,xlims=(xmin,xmax),xlabel="Time Step",
        yscale=:log10,ylims=(ymin,ymax),ylabel="Rt",
        title="Rt\n$(gname[1:end-5])")
    push!(result,f)

    #just looking at the deltas
    f = plot(leg=false)
    for k = 1:min(100,lastindex(inds)-1)
        #make deltas 
        delta_rt = rts[inds[k]:inds[k+1]-1]
        for i = lastindex(delta_rt):-1:2
            delta_rt[i] -= delta_rt[i-1]
        end
        plot!(f,ts[inds[k]:inds[k+1]-1],delta_rt,color=2,alpha=min(10/k,1),
            linestyle=:dash)
    end
    plot!(f,
        xscale=:log10,xlims=(xmin,105),xlabel="Time Step",
        ylabel="Change in Rt",
        title="Change in Rt\n$(gname[1:end-5])")
    abs_ylim = max(abs.(extrema(ylims(f)))...)
    plot!(f,
        ylims=(-abs_ylim,abs_ylim))
    push!(result,f)
    return result
end

function rt_plots1(gname::String,beta::Float64=0.05,qpercent::Int=0)
    
    if qpercent==0
        ts,rts,contagious_rts,exposed_rts = load_rtdata(gname,beta)
    else
        ts,rts,contagious_rts,exposed_rts = load_rtdata(gname,beta,qpercent)
    end
    inds = findall(ts.==1)

    result = []

    # making plots based on the batched data
        #contagious
    inds = findall(contagious_rts.>=0)
    batched_data = batchdata(ts[inds],contagious_rts[inds])
    inds = findall(length.(batched_data).>0)
    ys = mean.(batched_data[inds])

    f = scatter(inds,ys,xscale=:log10,xlims=(0.9,1e3),leg=false,
        title="Average Contagious Rt",
        markerstrokewidth=0)
    push!(result,f)

    #variance plot 
    ys = var.(batched_data[inds])
    f = scatter(inds,ys,xscale=:log10,xlims=(0.9,1e3),leg=false,
        title="Variance Contagious Rt",
        markerstrokewidth=0)
    push!(result,f)


        #exposed    
    inds = findall(exposed_rts.>=0)
    batched_data = batchdata(ts[inds],exposed_rts[inds])
    inds = findall(length.(batched_data).>0)
    ys = mean.(batched_data[inds])
    
    f = scatter(inds,ys,xscale=:log10,xlims=(0.9,1e3),leg=false,
        title="Average Exposed Rt",
        markerstrokewidth=0)
    push!(result,f)

    #variance plot 
    ys = var.(batched_data[inds])
    f = scatter(inds,ys,xscale=:log10,xlims=(0.9,1e3),leg=false,
        title="Variance Exposed Rt",
        markerstrokewidth=0)
    push!(result,f)

        #infected Rt
    inds = findall(rts.>=0)
    batched_data = batchdata(ts[inds],rts[inds])
    inds = findall(length.(batched_data).>0)
    ys = mean.(batched_data[inds])
    
    f = scatter(inds,ys,xscale=:log10,xlims=(0.9,1e3),leg=false,
        title="Average Infected Rt",
        markerstrokewidth=0)
    push!(result,f)    

    #variance plot 
    ys = var.(batched_data[inds])
    f = scatter(inds,ys,xscale=:log10,xlims=(0.9,1e3),leg=false,
        title="Variance Infected Rt",
        markerstrokewidth=0)
    push!(result,f)

    #heatmaps 
        #infected 
    inds = findall(rts.>0)
    f = myncpplot1(ts[inds],rts[inds],xlims=(0.9,1e3),ylims=(1e-5,1e3))
    plot!(f,title="Rt via all Infected nodes",xlabel="Time",ylabel="Rt")
    push!(result,f)

        #contagious
    inds = findall(contagious_rts.>0)
    f = myncpplot1(ts[inds],contagious_rts[inds],xlims=(0.9,1e3),ylims=(1e-5,1e3))
    plot!(f,title="Rt via all Contagious nodes",xlabel="Time",ylabel="Rt")
    push!(result,f)

        #exposed
    inds = findall(exposed_rts.>0)
    f = myncpplot1(ts[inds],exposed_rts[inds],xlims=(0.9,1e3),ylims=(1e-5,1e3))
    plot!(f,title="Rt via Exposed nodes",xlabel="Time",ylabel="Rt")
    push!(result,f)
    return result 
end

# gname = getgnames("flickr","input/graphs/")[1]
# gname = "rewired-10000.0-$gname"

# res = rt_plots1(gname,rt_param_dict[canonical_graph_name(gname)])


# res[1]
# res[2]
# res[3]

# res[4]
# res[5]
# res[6]

# res[7]
# res[8]
# res[9]





#using Makie to do this instead #needs Julia v1.8
using Makie
using CairoMakie
CairoMakie.activate!()

k = 3

delta_rts(ts,rts)
delta_rts(ts,rts)[inds[1]:inds[2]-1]
sliding_delta_rts = sliding_data_window(delta_rts(ts,rts)[inds[1]:inds[2]-1],20)
lastindex(sliding_delta_rts)

function rt_plots1(gname::String)

    result = []

    ts,rts = load_rtdata(gname)
    inds = findall(ts.==1)
    
    # plot type 1 #vanilla plot w/ alpha thresholding
    f = Figure()
    ax = Axis(f[1, 1],xlabel="Time Step",ylabel="Rt",
        title = "Rt\n$(gname[1:end-5])",
        xscale=log10)
        # yscale=log10)
    for k = 1:100
        tmprt = rts[inds[k]:inds[k+1]-1]
        Makie.lines!(ts[inds[k]:inds[k+1]-1],tmprt,linestyle=:dash,color=(:red,min(10/k,1)))
    end
    Makie.xlims!(ax,(0.9,700))
    # Makie.ylims!(ax,(1e-6,1e3))
    push!(result,f)

    #just looking at the deltas
    f = Figure()
    ax = Axis(f[1, 1],xlabel="Time Step",ylabel="Change in Rt",
        title = "Change in Rt\n$(gname[1:end-5])",
        xscale=log10,
        yscale=Makie.Symlog10(10.0),
        yticks = [-1000, -100, -10, 0, 10, 100, 1000])
    for k = 1:100
        delta_rt = rts[inds[k]:inds[k+1]-1]
        for i = lastindex(tmprt):-1:2
            delta_rt[i] -= delta_rt[i-1]
        end
        #aggregate over sliding window 
        final_rt = sliding_data_window(delta_rt,20)
        Makie.lines!(ts[inds[k]:inds[k+1]-1],final_rt,linestyle=:dash,color=(:red,min(10/k,1)))
    end
    Makie.xlims!(ax,(0.9,105))
    Makie.ylims!(ax,(-1e2,1e2))
    push!(result,f)
    return result
end



gname = gnames[27]
f,fdelta = rt_plots1(gname)
f
fdelta
#save figs 
figdst = "scratch/figures/rt/"
Makie.save(joinpath(figdst,"rt-1000-seir-0.05-0.05-$(gname[1:end-5]).png"),f)
Makie.save(joinpath(figdst,"delta-rt-1000-seir-0.05-0.05-$(gname[1:end-5]).png"),fdelta)

#=
# plot type 1 #vanilla plot w/ alpha thresholding
f = Plots.plot(leg=false)
for k = 1:100
    tmprt = rts[inds[k]:inds[k+1]-1]
    Plots.plot!(f,ts[inds[k]:inds[k+1]-1],tmprt,color=2,alpha=min(10/k,1),
        linestyle=:dash)
end
Plots.plot!(f,
    xscale=:log10,xlims=(xmin,xmax),xlabel="Time Step",
    yscale=:log10,ylims=(ymin,ymax),ylabel="Rt",
    title="Rt\n$(gname[1:end-5])")


#just looking at the deltas
f = Plots.plot(leg=false)
for k = 1:100
    tmprt = rts[inds[k]:inds[k+1]-1]
    for i = lastindex(tmprt):-1:2
        tmprt[i] -= tmprt[i-1]
    end
    Plots.plot!(f,ts[inds[k]:inds[k+1]-1],tmprt,color=2,alpha=min(10/k,1),
        linestyle=:dash)
end
plot!(f,
    xscale=:log10,xlims=(0.9,105),xlabel="Time Step",
    ylabel="Change in Rt",
    title="Change in Rt\n$(gname[1:end-5])")
=#
#batch the data together based on t
function batchdata(xs::Vector{T} where T,ys::Vector{S} where S) 
    batched_data = Vector{typeof(ys)}()
    for i=1:Int(maximum(xs))
        push!(batched_data,Vector{Float64}())
    end
    for (ind,val) in enumerate(xs)
        push!(batched_data[val],ys[ind])
    end
    return batched_data
end

function variance_plots(gname::String)
    f = Plots.plot()
    for g in ["er-10000.0-$gname" "rewired-10000.0-$gname" gname]
        ts,rts = load_rtdata(g)
        batched_rts = batchdata(ts,rts)

        vars = var.(batched_rts)
        ind = findfirst(isnan.(var.(batched_rts)))
        i2 = findfirst(vars.==0)
        if i2!=nothing || ind!=nothing 
            if ind != nothing && i2 != nothing 
                ind = min(i2,ind)
            elseif i2 != nothing 
                ind = i2
            end
        else 
            ind = lastindex(sliding_data)+1
        end

        if startswith(g,"er-")
            f = scatter(1:ind-1,vars[1:ind-1],label="ER",leg=:topright,
                markerstrokewidth=0,color=3)
        else 
            if startswith(g,"rewired")
                lab = "CL"
                color = 1
            else
                lab = "orig"
                color = 2
            end
            scatter!(f,1:ind-1,vars[1:ind-1],xscale=:log10,
                label=lab,color=color,#yscale=:log10,
                markerstrokewidth=0)    
        end
    end
    Plots.plot!(f,title=gname[1:end-5],xlabel="Time",ylabel="Variance in Rt")
    return f
end
variance_plots(gnames[1])
variance_plots(gnames[4])
variance_plots(gnames[7])
variance_plots(gnames[10])
variance_plots(gnames[13])
variance_plots(gnames[16])
variance_plots(gnames[19])
variance_plots(gnames[22])
variance_plots(gnames[25])
figdst = "scratch/figures/rt/"
for k=1:Int(lastindex(gnames)/3)
    gname = gnames[3*(k-1)+1]
    f = variance_plots(gname)
    Plots.savefig(f,joinpath(figdst,"rt-var-1000-seir-0.05-0.05-$(gname[1:end-5]).png"))
end

Makie.save(joinpath(figdst,"rt-1000-seir-0.05-0.05-$(gname[1:end-5]).png"),f)
Makie.save(joinpath(figdst,"delta-rt-1000-seir-0.05-0.05-$(gname[1:end-5]).png"),fdelta)


function delta_rts(ts::Vector{Int},rts::Vector{Float64})
    delta_rt = rts[:]
    for k = lastindex(rts):-1:2
        if ts[k]!=1
            delta_rt[k]-=delta_rt[k-1]
        end
    end
    return delta_rt
end

function sliding_data_window(batched_data::Vector{Vector{Float64}},windowsize::Int=10)
    mintime = windowsize
    result = Vector{Vector{Float64}}()
    for t = mintime:lastindex(batched_data)
        push!(result,vcat(batched_data[t-windowsize+1:t]...))
    end
    return result
end

function sliding_data_window(data::Vector{Float64},windowsize::Int=10)
    mintime = windowsize
    result = Vector{Vector{Float64}}()
    for t = mintime:lastindex(data)
        push!(result,vcat(data[t-windowsize+1:t]...))
    end
    return result
end

function variance_delta_plots(gname::String,windowsize::Int=10)
    f = Plots.plot()
    for g in ["er-10000.0-$gname" "rewired-10000.0-$gname" gname]
        ts,rts = load_rtdata(g)
        delta_rt = delta_rts(ts,rts)

        batched_delta_rt = batchdata(ts,delta_rt)
        
        sliding_data = sliding_data_window(batched_delta_rt,windowsize)
        
        vars = var.(sliding_data)
        ind = findfirst(isnan.(vars))
        i2 = findfirst(vars.==0)
        if i2!=nothing || ind!=nothing 
            if ind != nothing 
                ind = min(i2,ind)
            else
                ind = i2
            end
        else 
            ind = lastindex(sliding_data)+1
        end
        
        
        if startswith(g,"er-")
            f = scatter(1:ind-1,vars[1:ind-1],label="ER",leg=:topright,
                markerstrokewidth=0,color=3)
        else 
            if startswith(g,"rewired")
                lab = "CL"
                color = 1
            else
                lab = "orig"
                color = 2
            end
            scatter!(f,1:ind-1,vars[1:ind-1],xscale=:log10,
                label=lab,color=color,#yscale=:log10,
                markerstrokewidth=0)    
        end
    end
    Plots.plot!(f,title="$(gname[1:end-5])\nSliding Window Size: $windowsize",xlabel="Time",ylabel="Delta(Rt) Variance",
        ylims=(1e-15,ylims(f)[2]))
    return f
end
variance_delta_plots(gnames[1],20)


# figdst = "scratch/figures/rt/"
# for k=1:Int(lastindex(gnames)/3)
#     gname = gnames[3*(k-1)+1]
#     f = variance_delta_plots(gname,20)
#     Plots.savefig(f,joinpath(figdst,"rt-var-delta-20-1000-seir-0.05-0.05-$(gname[1:end-5]).png"))
# end

 

#Makie code 
f = Figure()
ax = Axis(f[1, 1],xlabel="Time Step",ylabel="Change in Rt",
    title = "Change in Rt\n$(gname[1:end-5])",
    xscale=log10,
    yscale=Makie.Symlog10(10.0))
for k = 1:100
    tmprt = rts[inds[k]:inds[k+1]-1]
    for i = lastindex(tmprt):-1:2
        tmprt[i] -= tmprt[i-1]
    end
    Makie.lines!(ts[inds[k]:inds[k+1]-1],tmprt,linestyle=:dash,color=(:red,min(10/k,1)))
end
Makie.xlims!(ax,(0.9,150))
f





#only look at increases in Rt
f = plot(leg=false)
for k = 1:100
    tmprt = rts[inds[k]:inds[k+1]-1]
    for i = lastindex(tmprt):-1:2
        tmprt[i] -= tmprt[i-1]
    end
    pos_deltas = findall(tmprt.>1e-1)
    data = tmprt[pos_deltas]
    #emphasize stronger flucuations
    scatter!(f,ts[inds[k]:inds[k+1]-1][pos_deltas],data,color=2,alpha=min(10/k,1),
        linestyle=:dash,markerstrokewidth=0)
end
plot!(f,xscale=:log10,xlims=(0.9,105))


=#