#script for producing figure comparing local structure (from ncp data) to error in uplots


# will make 3 plots.
# first explains measures for variance in local structure 
# second shows local structure versus error for CM
# second shows local structure versus error for GNP 
# one for CM and another for GNP 
mainDir = joinpath("/",split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))]...)
include(joinpath(mainDir,"code/graph-io.jl"))
include(joinpath(mainDir,"code/data-io.jl")) 

using DelimitedFiles
using DataFrames
using Measures
using Pkg
# using PerceptualColourMaps
using LaTeXStrings
using StatsBase
using ProgressMeter
using Plots
##
gr()
# ENV["GRDIR"] = ""
# Pkg.build("GR")


"""
    load_ncp(gname::String,datadir::String="pipeline/data/",normalizex::Bool=true;ncptype::String="epi")

loads in a single precomputed ncp. if normalizex is true, set size is made symmetric (min(size,nnode-size))
"""
function load_ncp(gname::String,datadir::String="pipeline/data/",normalizex::Bool=true;ncptype::String="epi")

  base_graph = canonical_graph_name(gname)
  if ncptype == "epi"
    fpath = joinpath(datadir,base_graph[1:end-5],"ncpdata","ncpinfo-epidemic-subsampling-4-seir-50000-$(gname[1:end-5]).txt")
  else
    fpath = joinpath(datadir,base_graph[1:end-5],"ncpdata","ncpinfo-$(gname[1:end-5]).txt")
  end

  ncp,headerinfo = readdlm(fpath,',',header=true)
  ncp = DataFrame(ncp,vec(headerinfo))
  
  
  if normalizex 
    n = graphSize(gname)
    x = map(x->min(ncp.size[x],n-ncp.size[x]),1:size(ncp,1))
  else 
    x = ncp.size
  end
  # x = ncp.size
  y = ncp.cond

  return x,y
end 


"""
    load_ncps(gnames::Vector,datadir::String="pipeline/data/",normalizex::Bool=true;ncptype::String="epi")

loads in ncps for every graph in gnames and extremal rewired variants. returns output as dict[gname]->(sizes,conds)
"""
function load_ncps(gnames::Vector,datadir::String="pipeline/data/",normalizex::Bool=true;ncptype::String="epi")
  data = Dict{String,Tuple}()
  @showprogress for gname in gnames 
    hs = ["rewired-10000.0-$gname" gname "er-10000.0-$gname"]
    for h in hs
      x,y = load_ncp(h,datadir,normalizex,ncptype=ncptype)

      data[h] = (x,y)
    end
  end
  return data
end

######### functions for computing "variance" of local structure 
#element wise variance function for nonzero values. 
function variance(data)
  return var(data[data.!=0])
end

#function for "area" of ncp
function count_nonzeros(X)
  c = 0
  for i in eachindex(X)
    c += Int(X[i]>0)
  end
  return c
end

#function for max-normalized nonzeros
function maxnormalize(X::Array)
  #max normalize columns 
  d = vec(mapslices(maximum,X;dims=1))
  Y = Float64.(X)
  for i in eachindex(d)
    if d[i]!=0
      Y[:,i]./=d[i]
    end
  end
  return Y
end

function count_nonzeros_threshold(X::Array,thresholds::Vector= [0.1;0.5;0.9])
  sort!(thresholds)
  @assert(all(0 .<= extrema(thresholds) .<= 1),"a least one threshold is not valid (<0 or >1)")
  X = maxnormalize(X)
  c = 0.0
  for i in eachindex(X)
    c+= sum(X[i] .< thresholds)
  end
  return c 
end


function variance_measure(f,gname::String;ncptype="epi")
    gnames = ["rewired-10000.0-$gname"; gname; "er-10000.0-$gname"]

    #predefine bin edges so we are talking about the same probability dist
    log_xedges = range(0,8,1601)  
    log_yedges = range(-6,0,1201)

    data = []
    for g in gnames
        x,y = load_ncp(g,ncptype=ncptype)

        h = StatsBase.fit(Histogram,(log10.(x),log10.(y)),(log_xedges,log_yedges))
        # @show size(h.weights)

        X = h.weights ./ sum(h.weights)
        push!(data,f(X))
    end
  
    return data
end

function variance_measure(f,ncpdata::Dict)
  #predefine bin edges so we are talking about the same probability dist
  log_xedges = range(0,8,1601)  
  log_yedges = range(-6,0,1201) 

  data = Dict{String,Float64}()
  for gname in keys(ncpdata)
      x,y = ncpdata[gname]
      h = StatsBase.fit(Histogram,(log10.(x),log10.(y)),(log_xedges,log_yedges))

      X = h.weights ./ sum(h.weights)
      data[gname] = f(X)
  end

  return data
end


#normalized variance measure 
function variance_measure1(f,ncpdata::Dict;
            sample_size_normalization::Bool=false)
  #predefine bin edges so we are talking about the same probability dist
  log_xedges = range(0,8,1601)  
  log_yedges = range(-6,0,1201) 

  data = Dict{String,Float64}()
  for gname in keys(ncpdata)
      x,y = ncpdata[gname]
      h = StatsBase.fit(Histogram,(log10.(x),log10.(y)),(log_xedges,log_yedges))

      X = h.weights ./ sum(h.weights)
      if sample_size_normalization #normalize by number of samples
        data[gname] = f(X)/lastindex(x)
      else
        data[gname] = f(X)
      end
  end

  return data
end

#option 2
#make coarser grid 
function variance_measure2(f,ncpdata::Dict;
            sample_size_normalization::Bool=false,
            nbins=(100,100))

  data = Dict{String,Float64}()
  for gname in keys(ncpdata)
      x,y = ncpdata[gname]
      h = StatsBase.fit(Histogram,(log10.(x),log10.(y)),nbins=nbins)

      X = h.weights ./ sum(h.weights)
      if sample_size_normalization #normalize by number of samples
        data[gname] = f(X)/lastindex(x)
      else
        data[gname] = f(X)
      end
  end

  return data
end

#coarser grid but same endpoints 
function variance_measure3(f,ncpdata::Dict;
            sample_size_normalization::Bool=false,
            nbins=(100,100))
  #predefine bin edges so we are talking about the same probability dist
  log_xedges = range(0,8,nbins[1])  
  log_yedges = range(-6,0,nbins[2]) 


  data = Dict{String,Float64}()
  for gname in keys(ncpdata)
      x,y = ncpdata[gname]
      h = StatsBase.fit(Histogram,(log10.(x),log10.(y)),(log_xedges,log_yedges))

      X = h.weights ./ sum(h.weights)
      if sample_size_normalization #normalize by number of samples
        data[gname] = f(X)/lastindex(x)
      else
        data[gname] = f(X)
      end
  end

  return data
end


#find expectation in y 
function variance_measure4(f,ncpdata::Dict;
            sample_size_normalization::Bool=false,
            nbins=(100,100))
  #predefine bin edges so we are talking about the same probability dist
  log_xedges = range(0,8,nbins[1])  
  log_yedges = range(-6,0,nbins[2]) 

  midpts = map(ind->(log_yedges[ind]+log_yedges[ind+1])/2,1:lastindex(log_yedges)-1)

  data = Dict{String,Float64}()
  for gname in keys(ncpdata)
      x,y = ncpdata[gname]           
      h = StatsBase.fit(Histogram,log10.(y),log_yedges)
      data[gname] = sum((h.weights ./sum(h.weights)).*midpts)
  end

  return data
end

#area to y=1
using CategoricalArrays

"""
    get_approximate_mins(x,y)

for data (x,y), segment the data into bins and get mins in each bin
"""
function get_approximate_mins(x,y)
  minx = Vector{Float64}()
  miny = Vector{Float64}()
  log10_xcuts = vcat(0,range(1,log10(maximum(x))+1e-5,501)) #1-10^6 uniformly on log scale 
  cv = CategoricalArrays.cut(log10.(x),log10_xcuts;allowempty=true)
  p = sortperm(cv)
  firstindex = 1
  while firstindex <= length(p)
    first = cv[p[firstindex]]
    lastindex = firstindex + 1
    while lastindex <= length(p) && cv[p[lastindex]] == first
      lastindex += 1
    end
    # get the index of the minimizing element of y
    imin = p[firstindex + argmin(@view y[p[firstindex:lastindex-1]]) - 1]
    #println(first, " ", firstindex, " ", lastindex, " ", imin)
    push!(minx, x[imin])
    push!(miny, y[imin])

    firstindex = lastindex # setup for next
  end
  return minx,miny
end

"""
    get_approximate_mins(x,y)

for data (x,y), segment the data into bins and get mins in each bin
"""
function get_approximate_mins(x,y,nbins=500)
  minx = Vector{Float64}()
  miny = Vector{Float64}()
  log10_xcuts = vcat(0,range(1,log10(maximum(x))+1e-5,nbins+1)) #1-10^6 uniformly on log scale 
  cv = CategoricalArrays.cut(log10.(x),log10_xcuts;allowempty=true)
  p = sortperm(cv)
  firstindex = 1
  while firstindex <= length(p)
    first = cv[p[firstindex]]
    lastindex = firstindex + 1
    while lastindex <= length(p) && cv[p[lastindex]] == first
      lastindex += 1
    end
    # get the index of the minimizing element of y
    imin = p[firstindex + argmin(@view y[p[firstindex:lastindex-1]]) - 1]
    #println(first, " ", firstindex, " ", lastindex, " ", imin)
    push!(minx, x[imin])
    push!(miny, y[imin])

    firstindex = lastindex # setup for next
  end
  return minx,miny
end


"""
    extraploate(minx::Vector,miny::Vector)

given minx and miny, extraploate out to nnodes along x using rightmost point
"""
function extrapolate!(minx::Vector,miny::Vector,nnodes::Union{Int,Nothing}=nothing)
  #get rightmost point
  xlast,ylast = minx[end],miny[end]

  if nnodes !== nothing && xlast<nnodes/2
    push!(minx,nnodes/2)
    push!(miny,ylast)
  end
  return true
end


"""
    approximate_auc(xs,ys)

approximate area under curve using trapezoidal rule
"""
function approximate_auc(xs,ys)
  @assert(lastindex(xs)==lastindex(ys),"inputs must have same length")
  
  p = sortperm(xs)

  area = 0
  for ind = 1:lastindex(xs)-1
    delta_x = xs[p[ind+1]]-xs[p[ind]]
    area += (0.5)*(delta_x)*(ys[p[ind]]+ys[p[ind+1]])
  end
  return area 
end


"""
    auc_measure(ncpdata::Dict;nbins=(61,81))

iterate over ncpdata and compute auc from y=1 to conductance values
"""
function auc_measure(ncpdata::Dict;graphSizes::Dict=Dict(),nbins=(61,81))
  data = Dict{String,Float64}()

  for gname in keys(ncpdata)
      x,y = ncpdata[gname]
     
      #get mins 
      minx,miny = get_approximate_mins(x,y)

      #extrapolate
      nnodes = graphSizes[canonical_graph_name(gname)]
      extrapolate!(minx,miny,nnodes)

      #compute auc from y=1 to conductance minimum
      auc = approximate_auc(minx, 1 .-miny) #auc from y=1 to min conductance
      data[gname] = auc/(nnodes/2) 
  end
  return data 
end


"""
    log_auc_measure(ncpdata::Dict;nbins=(61,81))

iterate over ncpdata and compute auc from y=1 to conductance values
"""
function log_auc_measure(ncpdata::Dict;graphSizes::Dict=Dict(),nbins=(61,81))
  data = Dict{String,Float64}()

  for gname in keys(ncpdata)
      x,y = ncpdata[gname]
      
      #get mins 
      minx,miny = get_approximate_mins(x,y)

      #extrapolate
      nnodes = graphSizes[canonical_graph_name(gname)]
      extrapolate!(minx,miny,nnodes)

      #compute auc from y=1 to conductance minimum
      # auc = approximate_auc(log10.(minx), log10.(1 .-miny)) #auc from y=1 to min conductance
      auc = approximate_auc(log10.(minx), -1 .*log10.(1 .-miny))
      data[gname] = auc/log10(nnodes/2) #normalizing 
  end
  return data 
end


# #testing 
# gname = gnames[2]
# x,y = epi_ncpdata[gname]
# #get mins 
# minx,miny = get_approximate_mins(x,y)
# minx[end]
# nnodes = nnodes_dict[gname]
# extrapolate!(minx,miny,nnodes)
# if minx[1]>1
#   minx = vcat(1.0,minx)
#   miny = vcat(1.0,miny)
# end

# minx./=(nnodes/2)

# f = plot(minx[1:end],miny[1:end],
#       xscale=:log10,
#       yscale=:log10,
#       leg=false,
#       xlabel="Fraction of Maximum NCP Size",
#       ylabel="Minimum Condutance",
#       c=:black,
#       dpi=1000,
#       ylims=(1e-4,1.0),
#       framestyle=:none,
#       linewidth=5,
#       background_color = :transparent)

# plot!(f,[minx[end-1:end]],[miny[end-1:end]],
#       linestyle=:dash,
#       markershape=:circle,
#       color=:red,
#       markerstrokewidth=0,
#       markersize=4)
# Plots.savefig(f,"/p/mnt/scratch/network-epi/code/paper-figs/example-figs/area-ncp-computation-pt2.png")
# plot!(f,minx,ones(lastindex(minx)),fillrange=miny.+3e-5,
#     c=:green, alpha=0.25)
# Plots.savefig(f,"/p/mnt/scratch/network-epi/code/paper-figs/example-figs/area-ncp-computation-pt3.png")
# # auc_measure(acl_ncpdata,graphSizes=nnodes_dict)[gname]

##bias towards sets of low conductance
#proportion of area below canonical_graph_name
function get_weight_matrix(gname::String,ncpdata::Dict,nbins=(121,161))
  x,y = ncpdata[gname]
  log_xedges = range(0,8,nbins[1])  
  log_yedges = range(-6,0,nbins[2]) 
  
  h = StatsBase.fit(Histogram,(log10.(x),log10.(y)),(log_xedges,log_yedges))
  X = h.weights ./ sum(h.weights)
  return X'
end


"""
    get_coverage(ncpdata::Dict,nbins=(1201,1601))

returns conductance needed to cover th of the probability ditribution along conductance
"""
function get_coverage(ncpdata::Dict;th=0.5,nbins=(1201,1601))
  
  @assert(0<=th<=1,"must have 0 <= th <= 1")

  log_xedges = range(0,8,nbins[1])  
  log_yedges = range(-6,0,nbins[2]) 

  data = Dict()
  for gname in keys(ncpdata)
    x,y = ncpdata[gname]
    
    h = StatsBase.fit(Histogram,(log10.(x),log10.(y)),(log_xedges,log_yedges))
    X = h.weights ./ sum(h.weights)
    X = X' #correct orientation 

    csum = cumsum(vec(sum(X;dims=2)))
    # ind = lastindex(csum)-findfirst(csum.>=th)
    ind = findfirst(csum.>=th)

    data[gname] = 10^(collect(log_yedges)[ind+1])
  end
  return data
end

## end variance functions 

#=
#testing  
log_xedges = range(0,8,1601)  
log_yedges = range(-6,0,1201) 

x,y = epi_ncpdata["modmexico-city.smat"]
h = StatsBase.fit(Histogram,(log10.(x),log10.(y)),(log_xedges,log_yedges))
X = h.weights ./ sum(h.weights)

Y = maxnormalize(X)
mapslices(extrema,X,dims=1)
mapslices(extrema,Y,dims=1)
=#

hs = ["commutes-all","mexico", "filtered", #row1
    "cn-moduillinois", "cn-Penn", "cn-modWiscon", #row2...
    "dblp","enron","anon", #row 3
    "cit-HepPh", "slashdot", "flickr", 
    "geometric","study-25-150","cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5.smat"]
    #  "rewired-10000.0-modmexico","rewired-10000.0-cn-moduillinois","er-10000.0-dblp-cc"]
gnames = map(h->h=="" ? "" : getgnames(h,"pipeline/graphs/")[1], hs)
titles = ["US Commutes", "Mexico City", "US Flows",
    "Illinois", "Penn", "Wisc.",    
    "Collaboration", "Email", "Anon",
    "Citation", "Slashdot","Flickr",    
    "Geometric", "Study", "LFR"] 

function load_diffusion_data(gname::String;
          betas::Vector{Float64}=vcat(1e-3:1e-3:1e-2,2e-2:1e-2:1e-1),
          method::String="seir")

  gnames = ["rewired-10000.0-$gname"; gname; "er-10000.0-$gname"]
  g = canonical_graph_name(gname)

  #load in data   
  data = Dict{String,Dict}() #gname -> beta -> total
  dloc = "pipeline/data/$(g[1:end-5])/diffusions/uniform/"
  for h in gnames
    h_data = Dict{Float64,Vector{Float64}}() 
    for beta in betas
      h_data[beta] = sum.(read_inf_data(h,dloc=dloc,beta=beta,dtype="tinfs",method=method))
    end
    data[h] = deepcopy(h_data)
  end

  return data
end


#### MAKE FINAL PLOT USING NNZ() AND FINALIZING
hs = ["commutes-all","mexico", "filtered", #row1
    "cn-moduillinois", "cn-Penn", "cn-modWiscon", #row2...
    "dblp","enron","anon", #row 3
    "cit-HepPh", "slashdot", "flickr", 
    "geometric","study-25-150","cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5.smat"]
    #  "rewired-10000.0-modmexico","rewired-10000.0-cn-moduillinois","er-10000.0-dblp-cc"]
gnames = map(h->h=="" ? "" : getgnames(h,"pipeline/graphs/")[1], hs)
titles = ["US-COMMUTES", "MEX", "US-FLO",
    "CN-UILL", "CN-PENN", "CN-WISC",    
    "COLLAB", "EMAIL", "FB-INT",
    "CIT", "SLA","FLICKR",    
    "GEO", "GEO-COMM", "LFR-RW"]


betas = [[0.003],[0.05],[0.03],
    [0.1],[0.1],[0.1],
    [0.1],[0.06],[0.1],
    [0.008],[0.04],[0.04],
    [0.06],[0.02],[0.175]]

# gnames = ["study-11-2023-0-noweak.smat",
# "study-11-2022-1.smat",
# "study-11-2023-1-longrange-1.smat",
# "study-11-2023-1-longrange-2.smat",
# "study-11-2023-1-longrange-3.smat",
# "study-11-2023-1-longrange-5.smat",
# "study-11-2023-1-longrange-8.smat",
# "study-11-2023-1-longrange-10.smat",
# "study-11-2023-1-longrange-12.smat",
# "study-11-2023-1-longrange-15.smat",
# "study-11-2022-10.smat",
# "study-11-2022-20.smat",
# "study-11-2022-30.smat",
# "study-11-2022-40.smat",
# "study-11-2022-50.smat"
# ]


# betas = [[0.03] for i =1:lastindex(gnames)]

# betas = [vcat(1e-3:1e-3:1e-2,2e-2:1e-2:1e-1) for i =1:lastindex(gnames)]

data = Dict{String,Dict}()
@showprogress for (ind,gname) in enumerate(gnames)
  tmp = load_diffusion_data(gname,betas = get_betas(gname))
  for k in keys(tmp)
    data[k] = tmp[k]
  end
end
#aggreate data 
aggregated_data = Dict{String,Dict}()
@showprogress for gname in gnames
  #load in graph once 
  nnodes = graphSize(gname)

  graph_variants = ["rewired-10000.0-$gname"; gname; "er-10000.0-$gname"]
  for h in graph_variants
    h_data = data[h]
    #process betas for a single graph 
    tmp = Dict()
    for beta in keys(h_data)
      tmp[beta] = aggregate_diffusion_data(hcat(h_data[beta]))./nnodes
    end
    #add to target 
    aggregated_data[h] = deepcopy(tmp)
  end
end 


# load ncp data
epi_ncpdata = load_ncps(gnames,ncptype="epi")
acl_ncpdata = load_ncps(gnames,ncptype="acl")

function qimpact_metric(data)
  return 1-( mean(data[2:end])/(data[1]) )
end

function qimpact_metric1(gnames::Vector,betas::Vector,diffusion_data::Dict)
  result = Dict{String,Float64}()
  for (i,gname) in enumerate(gnames)
    beta = betas[i][1]
    gdata = diffusion_data[gname][beta]

    result[gname] = 1-( mean(gdata[2:end])/(gdata[1]) )
    # result[gname] = 1-( gdata[end]/(gdata[1]) )
  end
  return result
end

function qimpact_metric2(gnames::Vector,betas::Vector,diffusion_data::Dict)
  result = Dict{String,Float64}()
  for (i,gname) in enumerate(gnames)
    beta = betas[i][1]
    gdata = diffusion_data[gname][beta]

    # result[gname] = 1-( mean(gdata[2:end])/(gdata[1]) )
    result[gname] = 1-( gdata[end]/(gdata[1]) )
  end
  return result
end


# #metric for impact of quarantine
# qimpact = qimpact_metric1(gnames,betas,aggregated_data)

# gname = getgnames("mex","input/graphs/")[1]
# tt = aggregated_data[gname][0.05]
# tt1 = aggregated_data["rewired-10000.0-$gname"][0.05]

# scatter(0:15,tt,label="orig")
# scatter!(0:15,tt1,label="CM")
# xlabel!("Quarantine Percent")
# ylabel!("Total Infections")

# 1-mean(tt[2:end])/tt[1]
# 1-mean(tt1[2:end])/tt1[1]

# 1-tt[end]/tt[1]
# 1-tt1[end]/tt1[1]


# 1-(mean(tt[2:end])/tt[1])
# 1-(mean(tt1[2:end])/tt1[1])
# qimpact[gname]


#area of local structure 
# epi_nnz_data = variance_measure3(count_nonzeros,epi_ncpdata,nbins=(80,200))
# acl_nnz_data = variance_measure2(count_nonzeros,acl_ncpdata,nbins=(80,200))

# variance_measure(count_nonzeros,acl_ncpdata)
# variance_measure1(count_nonzeros,acl_ncpdata,sample_size_normalization=true)
# variance_measure2(count_nonzeros,acl_ncpdata,sample_size_normalization=false,nbins=(100,200))



############ FINAL VERSION #############
#version 2 - use same network taxonomy but add annotations and use legend for taxonomy
# function base_plot_local_structure(gnames::Vector,xdata::Dict,ydata::Dict)
#   f = plot()
#   taxonomy = ["Human\nMobility","Sparsifed\nFriendship","Electronic\nInteractions","Distant\nNetworked\nInteractions","Synthetic"]

#   # cinds = round.(Int,range(1,lastindex(cgrad(:inferno)),5))
#   colors = [:green, :blue , :orange, :pink, :black]
#   markers = [:circle,:star5,:diamond, :utriangle, :star8]

#   # colors = cgrad(:inferno)[cinds]
#   for (i,gname) in enumerate(gnames)
#     cind = Int(ceil(i/3))
#     if i%3==1
#       scatter!(f,[xdata[gname]],[ydata[gname]],label=taxonomy[cind],
#         markerstrokewidth=0,markersize=8,color=colors[cind],
#         markershape=markers[cind])
#     else
#       scatter!(f,[xdata[gname]],[ydata[gname]],label=false,
#         markerstrokewidth=0,markersize=8,color=colors[cind],
#         markershape=markers[cind])
#     end
#   end
#   # plot!(f,legend=Symbol(:outer, :right),legendfontsize=7,
#   #       xlabel="Area of Epidemic NCP",
#   #       ylabel="Quarantine Impact on Total Infections")#,size=(800,500))
#   plot!(f,legend=:topleft,legendfontsize=7,
#         xlabel="Area of Epidemic NCP",
#         ylabel="Quarantine Impact on Total Infections")#,size=(800,500))
#   plot!(f,
#         top_margin=1Measures.mm,
#         right_margin=0Measures.mm)
#   return f
# end

# #previewing that we are doing the right thing
# f = base_plot_local_structure(gnames,acl_nnz_data,qimpact)
# plot!(f,xlabel="Normalized AUC of PPR NCP",xlims=(xlims(f)[1],1.05))

# f1 = base_plot_local_structure(gnames,epi_nnz_data,qimpact)
# xlabel!(f1,"Normalized AUC of Epidemic NCP")

# minx = xlims(f)[1]
# maxx = xlims(f)[2]
# minx = min(minx,xlims(f1)[1])
# maxx = max(maxx,xlims(f1)[2])
# plot(f,f1,layout=(2,1),size=(600,600),
#     xlims=(minx,maxx))


# #annotate plots and save 
# #acl/ppr data 
# newf = deepcopy(f)
# offsets = [
#   (-2e-2,2e-2),(1.5e-2,-4e-2),(0,-4e-2), #nMobility
#   (0,-4e-2)#=Uill=#,(0,1e-2)#=penn=#,(-4.5e-2,-1.5e-2)#=wisc=#, #sparsified
#   (2e-2,2e-2),(0,2e-2),(0,-4.5e-2), #Electronic
#   (0,-5e-2),(1e-2,2e-2),(0,3e-2), #Distant
#   (-3e-2,0),(-4e-2,1e-2),(2e-2,-4e-2), #Synthetic
# ]
# colors = [:green, :blue , :orange, :pink, :black]
# colors = [c for c in colors for i=1:3]

# for ind = 1:15
#   annotate!(newf,[([acl_nnz_data[gnames[ind]]+offsets[ind][1]],[qimpact[gnames[ind]]+offsets[ind][2]], (titles[ind], 8, colors[ind], :bottom))])
# end
# newf
# plot!(newf,dpi=1000)
# Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/example-figs/ppr-auc_local_structure-qimpact.png")



# #epidemic ncp data 
# newf = deepcopy(f1)
# offsets = [
#   (-2e-2,2e-2),(0,-4e-2),(0,-4e-2), #nMobility
#   (0,-4e-2),(-5e-2,1e-2),(6.5e-2,-1e-2), #sparsified
#   (-2e-2,2e-2),(0,2e-2),(0,-4.5e-2), #Electronic
#   (0,-5e-2),(1e-2,2e-2),(0,3e-2), #Distant
#   (-4e-2,0),(-3e-2,1e-2),(6.5e-2,-1.5e-2), #Synthetic
# ]
# colors = [:green, :blue , :orange, :pink, :black]
# colors = [c for c in colors for i=1:3]

# for ind = 1:15
#   annotate!(newf,[([epi_nnz_data[gnames[ind]]+offsets[ind][1]],[qimpact[gnames[ind]]+offsets[ind][2]], (titles[ind], 8, colors[ind], :bottom))])
# end
# plot!(newf,xlims=(minx,maxx))
# plot!(newf,dpi=1000)
# Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/example-figs/epi-auc_local_structure-qimpact.png")


####### ACL Version #######





#doing the fb ones out of curiousity
# fb_gnames = getgnames("cn","input/graphs/")
# # full_gnames = map(x->x[4:end],fb_gnames)
# # gnames = vcat(fb_gnames,full_gnames)
# filter!(x->!occursin("sbm",x),fb_gnames)

# fb_data = Dict{String,Dict}()
# @showprogress for gname in gnames
#   tmp = load_diffusion_data(gname,betas = get_betas(gname))
#   for k in keys(tmp)
#     fb_data[k] = tmp[k]
#   end
# end

# #aggreate data 
# fb_aggregated_data = Dict{String,Dict}()
# @showprogress for gname in gnames
#   #load in graph once 
#   nnodes = graphSize(gname)

#   graph_variants = ["rewired-10000.0-$gname"; gname; "er-10000.0-$gname"]
#   for h in graph_variants
#     h_data = data[h]
#     #process betas for a single graph 
#     tmp = Dict()
#     for beta in keys(h_data)
#       tmp[beta] = aggregate_diffusion_data(hcat(h_data[beta]))./nnodes
#     end
#     #add to target 
#     fb_aggregated_data[h] = deepcopy(tmp)
#   end
# end 

# # load ncp data
# fb_epi_ncpdata = load_ncps(fb_gnames,ncptype="epi")
# #metric for impact of quarantine
# fb_qimpact = qimpact_metric(fb_gnames,[[0.1] for i=1:lastindex(fb_gnames)],fb_aggregated_data)

# #area of local structure 
# fb_epi_nnz_data = variance_measure(count_nonzeros,fb_epi_ncpdata)

# f1 = plot()
# for (i,gname) in enumerate(fb_gnames)
#   if i==1
#     lab = "Sparsifed\nFriendship"
#   else
#     lab = false
#   end
#   scatter!(f1,[fb_epi_nnz_data[gname]],[fb_qimpact[gname]],
#     markerstrokewidth=0,markersize=8,c=:blue,marker=:star5,label=lab)
# end
# plot!(f1,
#     xlabel="Area of Epidemic NCP",
#     ylabel="Quarantine Impact on Total Infections")
    
# xticks!(f1,(1e4:1e4:7e4,["$(i)e4" for i=1:7]),size=(650,450))
# xlims!(f1,(8e3,3.2e4),leg=:right)


# #annotations
# println(fb_gnames)
# fb_titles = ["CN-UF","CN-UCF","CN-USF",
#             "CN-PENN","CN-FSU","CN-UNC",
#             "CN-ASTRO","CN-HARVARD","CN-STAN",
#             "CN-BERK", "CN-UILL","CN-WISC","CN-NOTRE-DAME"            
# ]

# offsets = [(0,-4e-2),(0,1e-2),(1e3,-3e-2),
#         (1e3,1e-2),(1.5e3,-2e-2),(0,-4e-2),
#         (0,1e-2),(0,1e-2),(0,1e-2),
#         (-1e3,-4e-2),(0,-3.5e-2),(-1.5e3,-2e-2),(0,-3.5e-2)]
# f2 = deepcopy(f1)

# for (i,gname) in enumerate(fb_gnames)
#   xloc = fb_epi_nnz_data[gname]+offsets[i][1]
#   yloc = fb_qimpact[gname]+offsets[i][2]
#   annotate!(f2,[([xloc],[yloc], (fb_titles[i], 8, :blue, :bottom))])
# end


# plot!(f2,dpi=1500)
# Plots.savefig("/p/mnt/scratch/network-epi/code/paper-figs/example-figs/local-structure-vs-qimpact-appendix.png")




########
# function base_plot_local_structure1(gnames::Vector,xdata::Dict,ydata::Dict)
#   f = plot()
#   taxonomy = ["Human\nMobility","Sparsifed\nFriendship","Electronic\nInteractions","Distant\nNetworked\nInteractions","Synthetic"]
#   taxonomy = vcat(taxonomy,["CM Rewired" for i=1:5])
#   # cinds = round.(Int,range(1,lastindex(cgrad(:inferno)),5))
#   colors = [:green, :blue , :orange, :pink, :black]
#   colors = vcat(colors,[:red for i=1:5])
#   markers = [:circle,:star5,:diamond, :utriangle, :star8]
#   markers = vcat(markers,[:pentagon for i=1:5])

#   # colors = cgrad(:inferno)[cinds]
#   for (i,gname) in enumerate(gnames)
#     cind = Int(ceil(i/3))
#     if i%3==1 && i<=18
#       scatter!(f,[xdata[gname]],[ydata[gname]],label=taxonomy[cind],
#         markerstrokewidth=0,markersize=8,color=colors[cind],
#         markershape=markers[cind])
#     else
#       scatter!(f,[xdata[gname]],[ydata[gname]],label=false,
#         markerstrokewidth=0,markersize=8,color=colors[cind],
#         markershape=markers[cind])
#     end
#   end
#   # plot!(f,legend=Symbol(:outer, :right),legendfontsize=7,
#   #       xlabel="Area of Epidemic NCP",
#   #       ylabel="Quarantine Impact on Total Infections")#,size=(800,500))
#   plot!(f,legend=:topleft,legendfontsize=7,
#         xlabel="Area of Epidemic NCP",
#         ylabel="Quarantine Impact on Total Infections")#,size=(800,500))
#   plot!(f,
#         top_margin=1Measures.mm,
#         right_margin=0Measures.mm)
#   return f
# end


#graph sizes and avgds
nnodes_dict = Dict()
avgd_data = Dict{String,Float64}()
deg_dist = Dict{String,Vector}()
@showprogress for gname in gnames
  @show(gname)
  A = loadGraph(gname,"input/graphs/")
  avgd_data[gname] = nnz(A)/size(A,1)
  nnodes_dict[gname] = size(A,1)
  deg_dist[gname] = vec(sum(A;dims=2))

  gname = "rewired-10000.0-$gname"
  @show(gname)
  A = loadGraph(gname,"pipeline/graphs/")
  avgd_data[gname] = nnz(A)/size(A,1)
  deg_dist[gname] = vec(sum(A;dims=2))
end

#eigvals
eig_data = readdlm("pipeline/data/dominant-eigvals-0.txt")
lam_data = Dict{String,Float64}()
for gname in gnames
  try
    ind = findfirst(eig_data[:,1].==gname)
    lam_data[gname] = eig_data[ind,2]

    gname = "rewired-10000.0-$gname"
    ind = findfirst(eig_data[:,1].==gname)
    lam_data[gname] = eig_data[ind,2]
  catch 
    println(gname)
  end
end

#raw data 
epi_nnz_data = auc_measure(epi_ncpdata,graphSizes=nnodes_dict)
acl_nnz_data = auc_measure(acl_ncpdata,graphSizes=nnodes_dict)

newgnames = vcat(gnames,map(x->"rewired-10000.0-$x",gnames))


# betas= [
#   0.03, #0
#   0.03, #1
#   2e-2, #L1
#   1e-2, #L2
#   9e-3, #L3
#   8e-3, #L5 
#   7e-3, #L8 
#   6e-3, #L10 
#   5e-3, #L12
#   4e-3, #L15
#   0.03, #10
#   0.03, #20 
#   0.03, #30 
#   0.03, #40 
#   0.03, #50
#   ]
# betas = vcat(betas,
#     [0.03,0.03,0.03,0.02,0.02,0.01,0.01,0.009,0.008,5e-3,
#     0.03,0.03,0.03,0.03,0.03])
newbetas = vcat(betas,betas)

#same color but different symbol 
function base_plot_local_structure2(gnames::Vector,xdata::Dict,ydata::Dict)
  f = plot()
  taxonomy = ["Human\nMobility","Sparsifed\nFriendship","Electronic\nInteractions","Distant\nNetworked\nInteractions","Synthetic"]
  taxonomy = vcat(taxonomy,["CM Rewired" for i=1:5])
  # cinds = round.(Int,range(1,lastindex(cgrad(:inferno)),5))
  cs = cgrad(:viridis)[Int.(range(1,236,6))]

  colors = reverse(cs[2:end])
  colors = vcat(colors,[cs[1] for i=1:5])
  # colors=[:cividis for i=1:10]

  markers = [:circle,:star5,:diamond, :utriangle, :star8]
  markers = vcat(markers,[:pentagon for i=1:5])

  # colors = cgrad(:inferno)[cinds]
  for (i,gname) in enumerate(gnames)
    cind = Int(ceil(i/3))
    if i%3==1 && i<=18
      scatter!(f,[xdata[gname]],[ydata[gname]],label=taxonomy[cind],
        markerstrokewidth=0,markersize=8,color=colors[cind],
        markershape=markers[cind])
    else
      scatter!(f,[xdata[gname]],[ydata[gname]],label=false,
        markerstrokewidth=0,markersize=8,color=colors[cind],
        markershape=markers[cind])
    end
  end
  # plot!(f,legend=Symbol(:outer, :right),legendfontsize=7,
  #       xlabel="Area of Epidemic NCP",
  #       ylabel="Quarantine Impact on Total Infections")#,size=(800,500))
  plot!(f,legend=:topleft,legendfontsize=7,
        xlabel="Area of Epidemic NCP",
        ylabel="Quarantine Impact on Total Infections")#,size=(800,500))
  plot!(f,
        top_margin=1Measures.mm,
        right_margin=0Measures.mm)
  return f
end

function log_auc_measure(ncpdata::Dict;graphSizes::Dict=Dict(),nbins=(61,81))
  data = Dict{String,Float64}()

  for gname in keys(ncpdata)
      x,y = ncpdata[gname]
      
      #get mins 
      minx,miny = get_approximate_mins(x,y)

      #extrapolate
      nnodes = graphSizes[canonical_graph_name(gname)]
      extrapolate!(minx,miny,nnodes)

      #compute auc from y=1 to conductance minimum
      auc = approximate_auc(log10.(minx), abs.(log10.(miny))) #auc from y=1 to min conductance
      data[gname] = auc/(log10(nnodes/2))
  end
  return data 
end

function base_plot_local_structure3(gnames::Vector,xdata::Dict,ydata::Dict)
  f = plot()
  colors = [:green, :blue , :orange, :pink, :black]
  colors = vcat(colors,[:red for i=1:5])

  markers = [:circle,:star5,:diamond, :utriangle, :star8]
  markers = vcat(markers,markers)#[:pentagon for i=1:5])

  # colors = cgrad(:inferno)[cinds]
  for (i,gname) in enumerate(gnames)
    cind = Int(ceil(i/3))
    if i%3==1 && i<=18
      scatter!(f,[xdata[gname]],[ydata[gname]],label=taxonomy[cind],
        markerstrokewidth=0,markersize=8,color=colors[cind],
        markershape=markers[cind])
    else
      scatter!(f,[xdata[gname]],[ydata[gname]],label=false,
        markerstrokewidth=0,markersize=8,color=colors[cind],
        markershape=markers[cind])
    end
  end
  # plot!(f,legend=Symbol(:outer, :right),legendfontsize=7,
  #       xlabel="Area of Epidemic NCP",
  #       ylabel="Quarantine Impact on Total Infections")#,size=(800,500))
  plot!(f,legend=:topleft,legendfontsize=7,
        xlabel="Area of Epidemic NCP",
        ylabel="Quarantine Impact on Total Infections")#,size=(800,500))
  plot!(f,
        top_margin=1Measures.mm,
        right_margin=0Measures.mm)
  return f
end

function entropy(data::Vector)
  #convert to probability distribution
  ps = counts(Int.(data))
  ps = ps./sum(ps) 
  
  #compute entropy
  result = 0
  for p in ps
    if p>0
      result += -p*log2(p)
    end
  end 
  return result/log2(lastindex(ps))
end

function kl_divergence(data::Vector)
  #convert to probability distribution
  ps = counts(Int.(data))
  ps = ps./sum(ps)   

  q = 1/lastindex(ps)

  result = 0
  for p in ps
    if p>0
      result += p*log2(p/q)
    end
  end 
  return result
end

epi0 = auc_measure(epi_ncpdata,graphSizes=nnodes_dict)
epi1 = log_auc_measure(epi_ncpdata,graphSizes=nnodes_dict)

qimpact = qimpact_metric1(newgnames,newbetas,aggregated_data)
ys = [qimpact[gname] for gname in newgnames]

x_auc = [epi1[gname] for gname in newgnames]
x_avgd = [avgd_data[gname] for gname in newgnames]
x_lam = [lam_data[gname] for gname in newgnames]


# scatter(x_auc, ys, color = :cividis, legend = false,
#   colorbar=true,markerstrokewidth=0,markersize=4,margins=5Measures.mm,
#   xscale=:log10)
# xlabel!("Local Structure (Normalized Epidemic NCP Area)",)
# ylabel!("Quarantine Impact")
# title!("Study Graphs - Area vs Impact")


# using PerceptualColourMaps


function strip_ytick_labels(f)
  y_ticks = yticks(f)[1]
  plot!(f,yticks=(y_ticks,["" for i=1:lastindex(y_ticks)]))
end


# Spearman rho stuff
## Centralize the codes for computing things...

using Distributions
function bootstrapped_mapped_correlation(x,y;ci=0.95,map=identity,ntrials=10^5)
  n = length(x)
  r = cor(map(x),map(y))
  ccvals = zeros(typeof(r), ntrials) 
  for i in 1:ntrials
    vals = rand(1:n, n)
    xr = map(x[vals])
    yr = map(y[vals])
    ccvals[i] = cor(xr,yr)
  end 
  cihalf = (1-ci)/2
  qvals = quantile(ccvals, (cihalf, 1-cihalf))
  return r, (qvals[1], qvals[2])
end 
function mapped_correlation_with_confidence(x,y;ci=0.95,map=identity)
  r = cor(map(x),map(y))
  zp = 0.5*(log(1+r) - log(1-r))
  qtile = (ci + 1)/2
  pval = quantile(Normal(), qtile)
  zpb = (zp-pval/sqrt(length(x)-3), zp+pval/sqrt(length(x)-3))
  rb = (exp.(2 .*zpb) .- 1)./(exp.(2 .* zpb) .+ 1)
  return r, rb
end 
function spearmanrho(x,y;ci=0.95,ntrials=10^5)
  map = x->invperm(sortperm(x))
  bsci = bootstrapped_mapped_correlation(x,y;map,ci,ntrials)[2] 
  r,aci = mapped_correlation_with_confidence(x,y;map,ci)
  return r, aci, bsci 
end 
function pearson_correlation(x,y;ci=0.95,ntrials=10^5)
  map = identity 
  bsci = bootstrapped_mapped_correlation(x,y;map,ci,ntrials)[2] 
  r,aci = mapped_correlation_with_confidence(x,y;map,ci)
  return r, aci, bsci 
end

spearmanrho(x_avgd,ys)
spearmanrho(x_lam,ys)
spearmanrho(x_auc,ys)

spearmanrho(data,ys)
# using Statistics
# using GLM

### explaination plot for appendix 
# gname = "modmexico-city.smat"
# x,y = epi_ncpdata[gname]
      
# #get mins 
# minx,miny = get_approximate_mins(x,y)
# minx = vcat(1,minx,minx[end])
# miny = vcat(1,miny,1)

# #extrapolate
# nnodes = nnodes_dict[gname]
# extrapolate!(minx,miny,nnodes)

# #compute auc from y=1 to conductance minimum
# nlines = 4
# ydim = (nlines+1)*35 + (nlines == 2 ? 1 : 0)

# f = plot(margin=-20Measures.mm,  topmargin=0Measures.mm, size=(250,ydim),
#               ylims=(10.0^(-nlines-1),2.0), framestyle=:none, xlims=(0.8,2e5), xscale=:log10, yscale=:log10)
# makegridlines!(f, [1,2,3,4,5], -(1:nlines),tickfontsize=15)
# plot!(minx,miny,xscale=:log10,yscale=:log10,leg=false,
#         linewidth=2,linecolor=:black)


# plot!(f,size=(600,400),dpi=1000,colorbar=true,
#       bottom_margin=-5Measures.mm)


# Plots.savefig(f,"/p/mnt/scratch/network-epi/code/paper-figs/example-figs/outline-mexicocity-ncp.png")
#         # fill=(1,:blue,0.3))
# # plot!(f,[nnodes/2;nnodes/2],[miny[end];ylims(f)[2]])

# f = plot(log10.(minx)./(log10(nnodes/2)),-log10.(miny),leg=false,
#           fill=(0,:black,0.3),
#           linewidth=2,linecolor=:black)

# # plot!(f,[log10(nnodes/2);log10(nnodes/2)],[ylims(f)[1];-log10(miny[end])])
# auc = approximate_auc(log10.(minx), abs.(log10.(miny))) #auc from y=1 to min conductance
# auc/(log10(nnodes/2))




epi1 = log_auc_measure(epi_ncpdata,graphSizes=nnodes_dict)
acl1 = log_auc_measure(acl_ncpdata,graphSizes=nnodes_dict)
qimpact = qimpact_metric1(newgnames,newbetas,aggregated_data)

tinfs = [aggregated_data[x[1]][x[2]][1] for x in zip(newgnames,newbetas)]
R0s = [lam_data[x[1]]*x[2]/5e-2 for x in zip(newgnames,newbetas)]

f = scatter(tinfs,leg=false,markerstrokewidth=0,markersize=5)
scatter!(twinx(),R0s,leg=false,c=:red,markerstrokewidth=0,markersize=3)




lam_data[newgnames[25]]*newbetas[25]/5e-2

x_avgd = [avgd_data[gname] for gname in newgnames]
x_lam = [lam_data[gname] for gname in newgnames]
x_epi = [epi1[gname] for gname in newgnames]
x_acl = [acl1[gname] for gname in newgnames]

ys = [qimpact[gname] for gname in newgnames]

f = plot()

scatter!(f,x_epi, ys, color = :cividis, legend = false,
  colorbar=true,markerstrokewidth=0,markersize=4,margins=5Measures.mm,
  xscale=:log10)
xlabel!("Local Structure (Normalized Epidemic NCP Area)",)
ylabel!(f,"Quarantine Impact")
title!(f,"Study Graphs - Area vs Impact")
#annote this 
gnames
ann = ["$x" for x in [0;1;"L1";"L2";"L3";"L5";"L8";"L10";"L12";"L15";10;20;30;40;50]]
ann = vcat(ann,["CM-$x" for x in ann])
gnames
ann_offsets = [
  (0,-5e-2), #0
  (0,-5e-2), #1
  (0,-5e-2), #longrange-1
  (0,-5e-2), #longrange-2
  (0,-5e-2), #longrange-3
  (0,-5e-2), #longrange-5
  (0,-5e-2), #longrange-8
  (0,-5e-2), #longrange-10
  (0,-5e-2), #longrange-12
  (0,-5e-2), #longrange-15
  (0,-5e-2), #10
  (0,-5e-2), #20
  (-2e-2,-4e-2), #30
  (0,-5e-2), #40
  (0,-5e-2), #50
  # (2e-2,5e-2), #cm-0
  # (2e-2,-4e-2), #cm-1
  # (6e-2,0), #cm-10
  # (6e-2,0), #cm-20
  # (6e-2,5e-2), #cm-30
  # (6e-2,2e-2), #cm-40
  # (6e-2,-1e-2), #cm-50
  ]

for (ind,_) in enumerate(gnames)
  xpos = x_epi[ind]+ann_offsets[ind][1]
  ypos = ys[ind]+ ann_offsets[ind][2]
  annotate!(f,xpos,ypos,text(ann[ind],7,:center))
end

plot!(f,xscale=:log10,ylims=(-0.05,1.1))
f



using DelimitedFiles
pdata = ["gname", "avgd", "lam1","normalized_epi_area","normalized_pagerank_area","qimpact"]
fname = "plotting_data-log_local_structure-qimpact-final-v1.txt"
open(joinpath(mainDir,fname), "w") do io
  data = vcat(reshape(pdata,(1,6)),hcat(newgnames, x_avgd, x_lam, x_epi, x_acl, ys ))
  writedlm(io, data)
end



### just using the raw data
## Test out the new plots for avgd & eigenvalue too.
# We want to see what plots like this look like
using StatsPlots, GLM
using DataFrames
using Measures
using DelimitedFiles 
data=readdlm(IOBuffer("""gname	avgd	lam1	normalized_epi_area	normalized_pagerank_area	qimpact
commutes-all.smat	100.36439790575916	261.73438567513426	0.833949362555367	0.9530841992573543	0.7547161080324665
modmexico-city.smat	21.197163536474328	459.4585956941147	1.458963479695048	1.8295552720656267	0.9228459405503989
covidflows-2020_08_31-filtered-20.smat	19.541350627379106	59.89837102287329	0.8171776811158877	1.2995608919290158	0.7269605217300695
cn-modUIllinois20.smat	2.994776602633283	28.50950762981281	0.5105698776973332	1.052354861172197	0.8909036387544063
cn-Penn94.smat	2.871672752288723	34.42067713491359	0.45177765501296374	1.1450561415234313	0.9113847079087134
cn-modWisconsin87.smat	2.6864014801110083	22.812197924567347	0.4486899360196691	0.9731211410332586	0.9084409544921495
dblp-cc.smat	6.328788541294007	75.12773836090547	0.7685988364986996	1.0001207956549498	0.5830648652511419
email-Enron.smat	10.731896961063628	118.41771488874616	0.4541478287917064	0.6609950109255293	0.386903021847418
anony-interactions-onemonthA-cc.smat	3.526345191040843	21.031051581299213	0.3126237624373247	0.8181940034421282	0.7674549466088011
cit-HepPh.smat	24.463474898985496	76.58116004020054	0.45785929808638887	0.5394930858563378	0.5671512714336859
Slashdot0811.smat	12.129782833505688	131.34180088134937	0.19874030474774024	0.26090292351762684	0.39566744888429317
flickr-links-sym.smat	19.048506084953033	1240.2994866017648	0.4491472369926089	0.8896117034306129	0.48610223432126065
geometric-100000-2d-1-20.smat	14.10362	17.656085282327656	1.2157989640131326	1.1049944836477057	0.9976684358096735
study-25-150.smat	31.91664	92.06024865135673	0.7295101154770703	0.6506595372753741	0.8541348350982128
cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5.smat	4.797183617495958	24.162599721644977	0.8723806399540358	1.28786653996601	0.660070091111056
rewired-10000.0-commutes-all.smat	100.12270597960871	175.25348068086936	0.054582274752836446	0.05175692844529757	0.05211631372895298
rewired-10000.0-modmexico-city.smat	20.93038131141462	255.7685088746265	0.10557350754680107	0.1039882202508409	0.2701438649614809
rewired-10000.0-covidflows-2020_08_31-filtered-20.smat	19.529410106023008	43.88795309441521	0.10840859203831028	0.08563910535908704	0.03612062272692662
rewired-10000.0-cn-modUIllinois20.smat	3.1752254390128143	23.449146613745263	0.272431428443058	0.43658152440850856	0.5819427887243289
rewired-10000.0-cn-Penn94.smat	3.073595887375543	27.914041547271573	0.22017877266109867	0.4652171643620516	0.6615423523447748
rewired-10000.0-cn-modWisconsin87.smat	2.8985094934751685	21.912171037203727	0.2566459896986516	0.4736445762270734	0.6318661595373567
rewired-10000.0-dblp-cc.smat	6.345307992608709	20.19870042821013	0.14539876612364522	0.23520622418713175	0.14129242386355867
rewired-10000.0-email-Enron.smat	10.594883362022912	89.05824863248971	0.11225245229807836	0.18286465701041194	0.3115520263969791
rewired-10000.0-anony-interactions-onemonthA-cc.smat	3.5571295891025163	14.470546019432424	0.15333221919008186	0.3992721865909504	0.7203684782140104
rewired-10000.0-cit-HepPh.smat	24.399697683206885	63.762784000300236	0.09705703901480553	0.1198527368811106	0.21201890853461935
rewired-10000.0-Slashdot0811.smat	12.030661840744571	117.60497207757928	0.10179864945286697	0.13099000357807059	0.4003668360736854
rewired-10000.0-flickr-links-sym.smat	18.835509083747723	733.3275303637621	0.08338540137219122	0.10326686317304863	0.5349237116977926
rewired-10000.0-geometric-100000-2d-1-20.smat	14.10216	15.071175571574479	0.08713105640078919	0.09133562859710925	0.002709943725980768
rewired-10000.0-study-25-150.smat	31.8948	46.59469195733523	0.06072250724950518	0.058139146603944834	0.00794432861171801
rewired-10000.0-cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5.smat	4.818749098997055	19.887516284706287	0.17196757970894258	0.3093574282527879	0.11154913561440294
"""), header=true)[1]

#ppr data and qimpact data 
xdata,ydata = data[:,5],data[:,end]

spearmanrho(data[:,2],data[:,end])
spearmanrho(data[:,3],data[:,end])
spearmanrho(data[:,4],data[:,end])
spearmanrho(data[:,5],data[:,end])

figs = []
fsize = (525,400)
ylab = "Quarantine Impact on Total Infections"

titles = ["Cmts","Mex","Flo",
        "UIll","Penn","Wisc",
        "Col","Ema","Int",
        "Cit","Sla","Fli",
        "G","GC","RWC"]
titles = vcat(titles,map(x->"CM-$x",titles))


#epidemic ncp data 
function final_area_impact_plot(area_data,impact_data)
  xs = log10.(Float64.(area_data))
  ys = impact_data
  mydata = DataFrame(map(v->(x=v[1],y=v[2]), zip(xs,ys)))

  r = lm( @formula(y ~ x), mydata)
  pred = predict(r, mydata, interval = :confidence, level = 0.995)
  #p = @df data scatter(:x, :y, leg = false)
  p = scatter(10 .^xs, ys,
        leg=false,
        size=fsize,
        color=:red3,
        xscale=:log10,
        markerstrokewidth=0,
        ylabel=ylab)

  # sort data on x
  xperm = sortperm(xs)
  xsort = xs[xperm]
  preds = pred[xperm,:]
  plot!(p, 10 .^xsort, preds.prediction, linewidth = 2, color=:lightgrey, z_order=:back,
          ribbon = (preds.prediction .- preds.lower, preds.upper .- preds.prediction))
  return p 
end

xdata,ydata = data[:,4],data[:,end]
p = final_area_impact_plot(xdata,ydata)

f = deepcopy(p)
ann_offsets = [
  (-1e-1,4e-2), (-1e-1,-4e-2), (-1e-1,-4e-2), #com, mex, flo,
  (0,-4e-2), (-8e-2,-3e-2), (0,4e-2), #uill, penn, wisc 
  (-1e-1,-4e-2), (0,-4e-2), (3e-2,0), #col, ema, int
  (2e-2,0), (2e-2,0), (-6e-2,0), #cit, sla, fli
  (-1e-2,4e-2), (0,4e-2), (0,-3e-2), #g, gc, rwc 

  #start of rewired variants 
  (-1e-3,4e-2), (-3.2e-2,1e-2), (-1e-2,4e-2), #com, mex, flo,
  (-5e-2,-4e-2), (-3e-2,4e-2), (2e-2,1e-2), #uill, penn, wisc 
  (0,4e-2), (0,4e-2), (0,4e-2), #col, ema, int
  (-1e-2,-3e-2), (-3e-2,0), (0,4e-2), #cit, sla, fli
  (-1e-2,-4e-2), (-1e-2,-4e-2), (0,-4e-2), #g, gc, rwc 
]

for ind in eachindex(data[:,1])
  x_off,y_off = ann_offsets[ind]
  annotate!(f,xdata[ind]+x_off,ydata[ind]+y_off,text(titles[ind],7,:left))
end

plot!(f,
    xlabel="Local Structure (Normalized Epidemic NCP Area)",
    margin=0Measures.mm,ylims=(-0.1,1.1))

push!(figs,f)


##average degree 
xdata,ydata = data[:,2],data[:,end]
p = final_area_impact_plot(xdata,ydata) 

f = deepcopy(p)
ann_offsets = [
  (-1e1,-4e-2), (-2,-4e-2), (-2,4e-2), #comm, mex, flo
  (0.2,-2e-2), (-1e-1,4e-2), (-3e-1,-3e-2), #uill, penn, wisc 
  (0.3,-3e-2), (-2,3e-2), (-4e-1,4e-2), #col, emi, int
  (0,4e-2), (1e0,-3e-2), (0,-4e-2), #cit, sla, fli
  (0,-4e-2), (0,4e-2), (8e-2,-3e-2), #g, gc, rwc 
  
  #start of rewired variants 
  (-2e1,4e-2), (0,4e-2), (0,4e-2), #comm, mex, flo
  (-2e-1,-4e-2), (-7e-1,3.5e-2), (2e-1,-1e-2), #uill, penn, wisc
  (0,4e-2), (0,-4e-2), (2e-1,1e-2),  #col, emi, int
  (0,-4e-2), (-5e-1,4e-2), (-3e0,4e-2), #cit, sla, fli
  (-3e0,4e-2), (1e0,-4e-2), (-5e-1,-4e-2), #g, gc, rwc
]


for ind in eachindex(data[:,1])
  x_off,y_off = ann_offsets[ind]
  annotate!(f,xdata[ind]+x_off,ydata[ind]+y_off,text(titles[ind],7,:left))
end

f

plot!(f,
    xlabel="Average Degree",
    margin=0Measures.mm,ylims=(-0.1,1.1))

push!(figs,f)


##lambda1 data 
xdata,ydata = data[:,3],data[:,end]
p = final_area_impact_plot(xdata,ydata) 

f = deepcopy(p)
ann_offsets = [
  (-1e1,-5e-2), (0,-5e-2), (-5,5e-2),  #comm, mex, flo
  (0.25,-4e-2), (0,4e-2), (-4,-3e-2), #uill, penn, wisc 
  (0.5,4e-2), (-10,-3e-2), (-1e0,4e-2), #col, emi, int
  (0,-4e-2), (1e1,-1e-2), (0,-4e-2), #cit, sla, fli
  (0,-5e-2), (0,5e-2), (-2,3e-2), #g, gc, rwc 
  
  #start of rewired variants 
  (-1e1,4e-2), (0,4e-2), (0,4e-2), #comm, mex, flo
  (-5,-4e-2), (0,-2e-2), (-9,-1e-2), #uill, penn, wisc
  (-5e0,4e-2), (-1e1,-4e-2), (0,-4e-2), # col, emi, int 
  (0,-4e-2), (-2e1,3e-2), (-1e0,4e-2), #cit, sla, fli
  (1e0,0), (3e0,0), (1e0,3e-2), # g, gc, rwc 
]

for ind in eachindex(data[:,1])
  x_off,y_off = ann_offsets[ind]
  annotate!(f,xdata[ind]+x_off,ydata[ind]+y_off,text(titles[ind],7,:left))
end

f

plot!(f,
    xlabel="Dominant Eigenvalue",
    margin=0Measures.mm,ylims=(-0.1,1.1))

push!(figs,f)


#laying things out 
newf = plot(figs...,layout=(1,3),size=(525*3,425),link=:y,
  ylims=(-0.1,1.05))

plot!(newf[1],
  left_margin=8Measures.mm,
  bottom_margin=9Measures.mm,
  right_margin=-5Measures.mm)

strip_ytick_labels(newf[2])
plot!(newf[2],ylabel="",
  right_margin=-5Measures.mm)
strip_ytick_labels(newf[3])
plot!(newf[3],ylabel="")

plot!(newf,dpi=1000)
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/example-figs/qimpact-v4.pdf")
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/example-figs/qimpact-v4.png")



#adding ppr ncp 
xdata,ydata = data[:,5],data[:,end]
p = final_area_impact_plot(xdata,ydata)

f = deepcopy(p)
ann_offsets = [
  (-8e-2,-4e-2), (-1e-1,-4e-2), (1e-1,0), #com, mex, flo,
  (-3e-2,-3e-2), (2e-2,3e-2), (-1e-1,3e-2), #uill, penn, wisc 
  (-1e-1,-4e-2), (0,-4e-2), (-5e-2,3e-2), #col, ema, int
  (-2e-2,-4e-2), (0,4e-2), (-6e-2,-4e-2), #cit, sla, fli
  (-1e-2,4e-2), (0,4e-2), (1e-1,0), #g, gc, rwc 

  #start of rewired variants 
  (-3e-3,4e-2), (-3.5e-2,0), (-1e-2,4e-2), #com, mex, flo,
  (-1e-1,-4e-2), (-3e-2,4e-2), (2e-2,-1e-2), #uill, penn, wisc 
  (-4e-2,4e-2), (-2e-2,4e-2), (-5e-2,4e-2), #col, ema, int
  (-1e-2,-3e-2), (-1e-2,4e-2), (0,4e-2), #cit, sla, fli
  (-1e-2,-4e-2), (-1e-2,-4e-2), (-2e-2,-4e-2), #g, gc, rwc 
]

for ind in eachindex(data[:,1])
  x_off,y_off = ann_offsets[ind]
  annotate!(f,xdata[ind]+x_off,ydata[ind]+y_off,text(titles[ind],7,:left))
end

plot!(f,
    xlabel="Normalized PPR NCP Area",
    margin=0Measures.mm,ylims=(-0.1,1.1))

push!(figs,f)


#layout first and last 
newf = plot(figs[1],figs[4],layout=(1,2),size=(525*2,425),link=:y,
  ylims=(-0.1,1.05))

plot!(newf[1],
  left_margin=8Measures.mm,
  bottom_margin=9Measures.mm,
  right_margin=-5Measures.mm)

strip_ytick_labels(newf[2])
plot!(newf[2],ylabel="",
  right_margin=-1Measures.mm)

plot!(newf[1],xlabel="Normalized Epidemic NCP Area")
plot!(newf[2],xlabel="Normalized PPR NCP Area")

plot!(newf,dpi=1000)
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/example-figs/qimpact-epi-vs-ppr.pdf")
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/example-figs/qimpact-epi-vs-ppr.png")







figs = []
xorig = Float64.(data[:,2])
mydata = DataFrame(map(v->(x=v[1],y=v[2]), zip(xorig,yorig)))

r = lm( @formula(y ~ x), mydata)
pred = predict(r, mydata, interval = :confidence, level = 0.995)

#p = @df data scatter(:x, :y, leg = false)
p = scatter(xorig, yorig, leg=false, size=fsize, color=:red3, markerstrokewidth=0,
  strokecolor=:white,
  markersize=5,
  xscale=:log10,
  xlabel="Average Degree",
  ylabel=ylab)

# sort data on x
xperm = sortperm(xorig)
xs = xorig[xperm]
preds = pred[xperm,:]
plot!(p, xs, preds.prediction, linewidth = 2, color=:lightgrey, z_order=:back,
        ribbon = (preds.prediction .- preds.lower, preds.upper .- preds.prediction))
push!(figs,p)





xorig = log10.(Float64.(data[:,2]))
mydata = DataFrame(map(v->(x=v[1],y=v[2]), zip(xorig,yorig)))

r = lm( @formula(y ~ x), mydata)
pred = predict(r, mydata, interval = :confidence, level = 0.995)

#p = @df data scatter(:x, :y, leg = false)
p = scatter(10 .^xorig, yorig, leg=false, size=fsize, color=:red3, markerstrokewidth=0,
  strokecolor=:white,
  markersize=5,
  xlabel="Average Degree",
  xscale=:log10,
  ylabel=ylab)

# sort data on x
xs = 10 .^xorig[xperm]
preds = pred[xperm,:]
plot!(p, xs, preds.prediction, linewidth = 2, color=:lightgrey, z_order=:back,
        ribbon = (preds.prediction .- preds.lower, preds.upper .- preds.prediction))
        
push!(figs,p)


plot(figs...,link=:y,size=(900,400),
  margins=5Measures.mm)



######################################## new approach ############
#normalize data so that data lies in [0,1] × [0,1]
# function normalize_ncpdata()
gnames[1]
gname = gnames[1]
x,y = epi_ncpdata[gname]

nnodes_dict[gname]

function test_area_function(gname::String,ncpdata::Dict,nnodes_dict::Dict)
  x,y = ncpdata[gname]
  
  #normalize data 
  x = log10.(x)./log10.((nnodes_dict[canonical_graph_name(gname)]/2))
  y = log10.(y)

  # fit histogram
  log_xedges = range(0,1,200)  
  log_yedges = range(-6,0,100) 

  h = StatsBase.fit(Histogram,(x,y),(log_xedges,log_yedges))
  X = h.weights ./ sum(h.weights)

  #expectation in y 
  f = heatmap(maxnormalize(Array(X')))
  ymidpts = map(x->(log_yedges[x]+log_yedges[x+1])/2,1:(lastindex(log_yedges)-1))
  Dy = diagm(0=>-1 .*ymidpts)
  sum(X)
  sum(X*Dy)
  return f,sum(X*Dy)
end

#normalize x 
x = log10.(x)./log10.((nnodes_dict[gname]/2))
y = log10.(y)

# fit histogram
log_xedges = range(0,1,200)  
log_yedges = range(-6,0,100) 

h = StatsBase.fit(Histogram,(x,y),(log_xedges,log_yedges))
X = h.weights ./ sum(h.weights)

#expectation in y 
heatmap(maxnormalize(Array(X')))
ymidpts = map(x->(log_yedges[x]+log_yedges[x+1])/2,1:(lastindex(log_yedges)-1))
Dy = diagm(0=>-1 .*ymidpts)
sum(X)
sum(X*Dy)


gname = gnames[12]
gname = "rewired-10000.0-$(gnames[12])"
f,metric = test_area_function(gname,epi_ncpdata,nnodes_dict)
metric
f





#=
study_data = readdlm(IOBuffer("""gname	avgd	lam1	normalized_epi_area	normalized_pagerank_area	qimpact
study-11-2023-0-noweak.smat	10.70584	23.30000354864058	1.2169331402638752	1.1264883987689756	0.9967272201700336
study-11-2022-1.smat	11.20332	23.31313096515302	0.5631727158998859	0.6926873894434257	0.9959306651540681
study-11-2023-1-longrange-1.smat	11.70272	23.336993803458835	0.49961042300126707	0.5285696389709046	0.9952144646042713
study-11-2023-1-longrange-2.smat	12.70684	23.39918892346967	0.3745440824248874	0.39879259679383217	0.9811816302820834
study-11-2023-1-longrange-3.smat	13.7008	23.423497009889264	0.3062758904988687	0.32151526672275854	0.9456899854252725
study-11-2022-10.smat	15.84708	27.63348194834611	0.6090788889448927	0.7095449392825621	0.9946336087358602
study-11-2022-20.smat	18.68984	34.71380978882812	0.7092518982906814	0.7737849821106941	0.9789576580500452
study-11-2022-30.smat	21.3154	44.33256336899045	0.9156362513695986	0.9073509342665756	0.9675760668649039
study-11-2022-40.smat	23.40704281712685	51.18788682330759	0.9415806956328907	0.932850080487554	0.9459398988568847
study-11-2022-50.smat	24.865809007422822	52.135264130556386	1.0586982295090233	0.9354061261486514	0.942086236239391
rewired-10000.0-study-11-2023-0-noweak.smat	10.70296	18.193779760106395	0.10588882261214484	0.11437836738017645	0.5085654845766922
rewired-10000.0-study-11-2022-1.smat	11.20036	18.45968831235343	0.10847040694566541	0.1133246902239799	0.492428889389842
rewired-10000.0-study-11-2023-1-longrange-1.smat	11.69988	18.774294529263443	0.09905331582698232	0.10948789350558998	0.4293016324210185
rewired-10000.0-study-11-2023-1-longrange-2.smat	12.70392	19.1133659863091	0.08641560440803832	0.09868855130712231	0.4172102762895017
rewired-10000.0-study-11-2023-1-longrange-3.smat	13.6974	19.615249650058253	0.07967143486685015	0.0967926890784145	0.39406485582623096
rewired-10000.0-study-11-2022-10.smat	15.84364	17.981711533549127	0.08471051402969175	0.08724634222510554	0.3556028689623606
rewired-10000.0-study-11-2022-20.smat	18.68468	21.584024267034557	0.07683291592229356	0.07474022996996421	0.16057304353645063
rewired-10000.0-study-11-2022-30.smat	21.30852	25.147729850885906	0.07373433515471392	0.07315456849421556	0.005091755054808522
rewired-10000.0-study-11-2022-40.smat	23.39779911964786	28.121165258676534	0.06808011860633857	0.06982669190483606	0.0028939828577616566
rewired-10000.0-study-11-2022-50.smat	24.856605510093836	30.14032791368516	0.06902546103926842	0.06864672962340532	0.0022449489649605248
"""))


data
study_data
data = vcat(data,study_data[2:end,:])

data

fsize = (525,400)
ylab = "Quarantine Impact on Total Infections"

titles = ["Cmts","Mex","Flo",
        "UIll","Penn","Wisc",
        "Col","Ema","Int",
        "Cit","Sla","Fli",
        "G","GC","RWC",
        "GC-0","GC-1","GC-10",
        "GC-20","GC-30","GC-40","GC-50"]
titles = vcat(titles,map(x->"CM-$x",titles))

yorig = Float64.(data[:,end])
#epidemic ncp data 
xorig = log10.(Float64.(data[:,4]))
mydata = DataFrame(map(v->(x=v[1],y=v[2]), zip(xorig,yorig)))

r = lm( @formula(y ~ x), mydata)
pred = predict(r, mydata, interval = :confidence, level = 0.995)
#p = @df data scatter(:x, :y, leg = false)
p = scatter(10 .^xorig, yorig,
      leg=false,
      size=fsize,
      color=:red3,
      xscale=:log10,
      markerstrokewidth=0,
      ylabel=ylab)


# sort data on x
xperm = sortperm(xorig)
xs = xorig[xperm]
preds = pred[xperm,:]
plot!(p, 10 .^xs, preds.prediction, linewidth = 2, color=:lightgrey, z_order=:back,
        ribbon = (preds.prediction .- preds.lower, preds.upper .- preds.prediction))

plot!(p,xlabel="Local Structure (Normalized Epidemic NCP Area)")

f = deepcopy(p)
ann_offsets = [
  (-1e-1,4e-2), (-1e-1,-4e-2), (-1e-1,-4e-2), #com, mex, flo,
  (0,-4e-2), (-8e-2,-3e-2), (0,4e-2), #uill, penn, wisc 
  (-1e-1,-4e-2), (0,-4e-2), (3e-2,0), #col, ema, int
  (2e-2,0), (2e-2,0), (-6e-2,0), #cit, sla, fli
  (-1e-2,4e-2), (-1e-1,-4e-2), (0,3e-2), #g, gc, rwc 

  #start of rewired variants 
  (-1e-3,4e-2), (-3.2e-2,1e-2), (-1e-2,4e-2), #com, mex, flo,
  (-5e-2,-4e-2), (-3e-2,4e-2), (2e-2,1e-2), #uill, penn, wisc 
  (0,4e-2), (0,4e-2), (0,4e-2), #col, ema, int
  (-1e-2,-3e-2), (-3e-2,0), (0,4e-2), #cit, sla, fli
  (-1e-2,-4e-2), (-1e-2,-4e-2), (0,-4e-2), #g, gc, rwc 
]
=#