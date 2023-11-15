# study testing
mainDir = "/p/mnt/scratch/network-epi/code/"
include(joinpath(mainDir,"graph-io.jl"))
include(joinpath(mainDir,"data-io.jl"))
include(joinpath(mainDir,"ncp/ncp-acl.jl"))

#code for network generation 
using NearestNeighbors, Distributions, SparseArrays, LinearAlgebra
using MatrixNetworks, Arpack
function spatial_graph_edges(n::Integer,d::Integer;degreedist=LogNormal(log(4),1))
  xy = rand(d,n)
  T = BallTree(xy)
  # form the edges for sparse
  ei = Int[]
  ej = Int[]
  for i=1:n
    deg = min(ceil(Int,rand(degreedist)),n-1)
    idxs, dists = knn(T, xy[:,i], deg+1)
    for j in idxs
      if i != j
        push!(ei,i)
        push!(ej,j)
      end
    end
  end
  return xy, ei, ej
end
"""
n nodes in dimension d (random [0,1] box)
degreedist = degree distribution function
p = number of expected edges per node to add at ranodm (default = 0 )
"""
function spatial_network(n::Integer, d::Integer; degreedist=LogNormal(log(3),1),
    p::Real=0.0)
  xy, ei, ej = spatial_graph_edges(n, d;degreedist=degreedist)
  A = sparse(ei,ej,1,n,n)
  C = A
  if p > 0
    C .+= triu(sprand(Bool, n,n,p/n),1)
    fill(C.nzval,1)
  end
  return max.(C,C'), xy
end

function adjust_graph(A::SparseMatrixCSC, xy;
    factor::Real = 0.5, longedges::Real = 0, tdegs::Vector )
  n = size(A,1)
  degs = vec(sum(A; dims=2))
  newxy = copy(xy)
  rowval = rowvals(A)
  ei = Int[]
  ej = Int[]
  # compute average length of edges to find long ones...
  threshlen = Inf
  if longedges > 0
    avglen = 0.0
    for i=1:size(A,1)
      for nzi = nzrange(A,i)
        j = rowval[nzi]
        avglen += sqrt((xy[1,i]-xy[1,j])^2 + (xy[2,i]-xy[2,j])^2) / nnz(A)
      end
    end
    threshlen = longedges*avglen
  end
  rdegs = copy(degs) # residual degrees
  for i=1:size(A,1)
    bigj = i
    myd = degs[i]
    for nzi = nzrange(A,i)
      j = rowval[nzi]
      if degs[j] > degs[bigj]
        bigj = j
      end
      if longedges > 0 && sqrt((xy[1,i]-xy[1,j])^2 + (xy[2,i]-xy[2,j])^2) > threshlen
        push!(ei, i)
        push!(ej, j)
        rdegs[i] -= 1
        if rand() < (degs[j]/(degs[i]+degs[j]))
          rdegs[j] -= 1
        end
      end
    end
    if bigj != i # then we have a bigger neighbor
      newxy[:,i] = factor*newxy[:,i] + (1-factor)*newxy[:,bigj]  + 0.001*randn(2)
    end
  end
  #@show minimum(rdegs)
  # find new edges
  # form the edges for sparse
  T = BallTree(newxy)
  for i=1:n
    #deg = ceil(Int, degs[i])
    deg = ceil(Int, 0.85*max(ceil(Int,rdegs[i]),0))
    deg += ceil(Int, max(0, 0.25*(tdegs[i]-degs[i])))
    #deg = degs[i]
    idxs, dists = knn(T, newxy[:,i], deg+1)
    for j in idxs
      if i != j
        push!(ei,i)
        push!(ej,j)
      end
    end
  end
  B = sparse(ei,ej,1,n,n)
  fill!(B.nzval, 1)
  #= fixup the nodes that are now disconnected...
  for these, we add an edge to a random node in the
  big cc based on degree
  connect. =#
  B = max.(B,B')
  cc = scomponents(B)
  newdegs = vec(sum(B;dims=2))
  bigcc = argmax(cc.sizes)
  for i=1:size(B,1)
    if cc.map[i] == bigcc
      #continue # don't do anything
      if tdegs[i] > newdegs[i]
        j = rand(1:size(B,1))
        if tdegs[j] > newdegs[j]
          push!(ei, i)
          push!(ej, j) # add some random edges...
        end
      end
    else
      # not in big CC, need to add an edge
      if rand() <= 2/cc.sizes[cc.map[i]] # add two edges from components...
        j = rand(1:size(B,1))
        while cc.map[j] != bigcc # while we aren't in the big cc...
          j = rand(1:size(B,1))
        end
        push!(ei,i)
        push!(ej,j)
      end
    end
  end
  D = sparse(ei,ej,1,n,n)
  fill!(D.nzval, 1)
  D = max.(D,D')
  #@assert is_connected(D)
  return D, newxy
end

function adjust_graph_new(A::SparseMatrixCSC, xy;
    factor::Real = 0.5, longedges::Real = 0, tdegs::Vector )
  n = size(A,1)
  degs = vec(sum(A; dims=2))
  newxy = copy(xy)
  rowval = rowvals(A)
  ei = Int[]
  ej = Int[]
  # compute average length of edges to find long ones...
  threshlen = Inf
  if longedges > 0
    avglen = 0.0
    for i=1:size(A,1)
      for nzi = nzrange(A,i)
        j = rowval[nzi]
        avglen += sqrt((xy[1,i]-xy[1,j])^2 + (xy[2,i]-xy[2,j])^2) / nnz(A)
      end
    end
    threshlen = longedges*avglen
  end
  rdegs = copy(degs) # residual degrees
  for i=1:size(A,1)
    bigj = i
    myd = degs[i]
    for nzi = nzrange(A,i)
      j = rowval[nzi]
      if degs[j] > degs[bigj]
        bigj = j
      end
      if longedges > 0 && sqrt((xy[1,i]-xy[1,j])^2 + (xy[2,i]-xy[2,j])^2) > threshlen
        push!(ei, i)
        push!(ej, j)
        rdegs[i] -= 1
        # if rand() < (degs[j]/(degs[i]+degs[j]))
        #   rdegs[j] -= 1
        # end
      end
    end
    if bigj != i # then we have a bigger neighbor
      newxy[:,i] = factor*newxy[:,i] + (1-factor)*newxy[:,bigj]  + 0.001*randn(2)
    end
  end
  #@show minimum(rdegs)
  # find new edges
  # form the edges for sparse
  T = BallTree(newxy)
  for i=1:n
    #deg = ceil(Int, degs[i])
    deg = ceil(Int, 0.85*max(ceil(Int,rdegs[i]),0))
    deg += ceil(Int, max(0, 0.15*(tdegs[i]-degs[i])))
    #deg = degs[i]
    idxs, dists = knn(T, newxy[:,i], deg+1)
    for j in idxs
      if i != j
        push!(ei,i)
        push!(ej,j)
      end
    end
  end
  B = sparse(ei,ej,1,n,n)
  fill!(B.nzval, 1)
  #= fixup the nodes that are now disconnected...
  for these, we add an edge to a random node in the
  big cc based on degree
  connect. =#
  B = max.(B,B')
  cc = scomponents(B)
  newdegs = vec(sum(B;dims=2))
  bigcc = argmax(cc.sizes)
  for i=1:size(B,1)
    if cc.map[i] == bigcc
      #continue # don't do anything
      if tdegs[i] > newdegs[i]
        j = rand(1:size(B,1))
        if tdegs[j] > newdegs[j]
          push!(ei, i)
          push!(ej, j) # add some random edges...
        end
      end
    else
      # not in big CC, need to add an edge
      if rand() <= 2/cc.sizes[cc.map[i]] # add two edges from components...
        j = rand(1:size(B,1))
        while cc.map[j] != bigcc # while we aren't in the big cc...
          j = rand(1:size(B,1))
        end
        push!(ei,i)
        push!(ej,j)
      end
    end
  end
  D = sparse(ei,ej,1,n,n)
  fill!(D.nzval, 1)
  D = max.(D,D')
  #@assert is_connected(D)
  return D, newxy
end

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

#code for probing those networks 
include(joinpath(mainDir,"fast-diffusion.jl"))

function histogram_to_data(h)
        #use midpts for hexbinning
        xres,yres = Vector{Float64}(),Vector{Float64}()

        xmids = map(a->(h.edges[1][a]+h.edges[1][a+1])/2,1:lastindex(h.edges[1])-1)
        ymids = map(a->(h.edges[2][a]+h.edges[2][a+1])/2,1:lastindex(h.edges[2])-1)

        for xind = 1:lastindex(xmids)
                for yind = 1:lastindex(ymids)   
                        wind = h.weights[xind,yind]
                        if wind>0
                                append!(xres,repeat([xmids[xind]],wind))
                                append!(yres,repeat([ymids[yind]],wind))
                        end
                end
        end
        return xres,yres 
end

function epi_sample(E::Union{EventEpidemic.SEIRData,EventEpidemic.SIRData},
                    ktrials::Int)
  
        #run ntrials and if we have more than kfailures, pick a different seed node and try again
        seed_node = rand(1:length(E.snodes))
        z = zeros(Float64,length(E.snodes))
        for trial_num = 1:ktrials
                l,E = EventEpidemic.epidemic(E,seed_node)
                netinfs, newinfs = EventEpidemic.get_infdata(l, E)
                c = 0
                while  sum(newinfs) < min(500,0.01*lastindex(z)) && c<5
                        c+=1  
                        l,E = EventEpidemic.epidemic(E,seed_node)
                        netinfs, newinfs = EventEpidemic.get_infdata(l, E)
                end
                #define weight fcn for nodes 
                weight_fcn(x) = E.snodes[x] ? 0.0 : 1/(1+E.itime[x])
                z.+=min.(E.itime,l+1)
        end
        return -z
end


function epi_sample_figure(A::SparseMatrixCSC,beta::Float64=1e-1)
    E = EventEpidemic.SEIRData(A,beta=beta)

    xedges = vcat(0:1:99,100:10:1000-1,1000:50:1e4-1,1e4:5e2:1e5)
    yedges = 10.0 .^(range(-6,0,120))
    dvec = vec(sum(A;dims=2))

    #initialize histogram 
    h = fit(Histogram,([0],[0]),(xedges,yedges))

    @showprogress for i=1:100
        xrank = epi_sample(E,3)
        p = MatrixNetworks.sweepcut(A,xrank);
        StatsBase.merge!(h,fit(Histogram,(1:lastindex(p.conductance),p.conductance),(xedges,yedges)))
    end

    #get data for hexbin plotting 
    xres,yres = histogram_to_data(h)

    xbounds,ybounds = (1.0,2e5),(2e-6,1.15)
    c = cgrad(cgrad(:viridis)[1:250])

    f = plot()
    # myhexbin_conditional!(f,xres,yres,xlims=xbounds,ylims=(1e-2,2),color=c,nbins=80)
    myhexbin2!(f,xres,yres,xlims=xbounds,ylims=(1e-3,2),color=c,nbins=80,normalize_color=true)
    myhexbin2!(f,xres,yres,xlims=xbounds,ylims=(1e-3,2),color=c,nbins=80,normalize_color=true)
    plot!(f,xticks=[1,10,100,1000,1e4,1e5],
        yticks=[1e-5,1e-4,1e-3,1e-2,1e-1,1e0])
    return f,xres,yres 
end

function ppr_sample_figure(A::SparseMatrixCSC,ntrials::Int=100)
    xedges = vcat(0:1:99,100:10:1000-1,1000:50:1e4-1,1e4:5e2:1e5)
    yedges = 10.0 .^(range(-6,0,120))
    dvec = vec(sum(A;dims=2))
    
    nnodes = lastindex(A,1)

    #initialize histogram 
    h = fit(Histogram,([0],[0]),(xedges,yedges))

    @showprogress for i=1:ntrials
        ppr_x = MatrixNetworks.personalized_pagerank(A,0.995,rand(1:lastindex(A,1)),1e-2);
        p = MatrixNetworks.sweepcut(A,ppr_x./dvec);
        StatsBase.merge!(h,fit(Histogram,(1:lastindex(p.conductance),p.conductance),(xedges,yedges)))
        xs = 1:lastindex(p.conductance)
        xs = map(x->min(x,nnodes-x),xs)
        StatsBase.merge!(h,fit(Histogram,(xs,p.conductance),(xedges,yedges)))
    end

    #get data for hexbin plotting 
    xres,yres = histogram_to_data(h)

    xbounds,ybounds = (1.0,2e5),(2e-6,1.15)
    c = cgrad(cgrad(:viridis)[1:250])

    f = plot()
    # myhexbin_conditional!(f,xres,yres,xlims=xbounds,ylims=(1e-2,2),color=c,nbins=80)
    myhexbin2!(f,xres,yres,xlims=xbounds,ylims=(1e-3,2),color=c,nbins=80,normalize_color=true)
    myhexbin2!(f,xres,yres,xlims=xbounds,ylims=(1e-3,2),color=c,nbins=80,normalize_color=true)
    plot!(f,xticks=[1,10,100,1000,1e4,1e5],
        yticks=[1e-5,1e-4,1e-3,1e-2,1e-1,1e0])
    return f,xres,yres 
end


# using Graphs

# nnodes = 50000
# ngroups = 50
# avgd = 20

# pdensity = zeros(ngroups,ngroups)
# for i=1:ngroups
#   cin = avgd/2+rand()*(avgd/2)
#   cout = avgd-cin
#   cout,cin = extrema([cin,cout])

#   cout = 0.05*cout

#   #diag
#   pdensity[i,i] = cin

#   #split outdegree among all groups 
#   v = rand(ngroups-1)
#   v./=sum(v)

#   v = cout*v
#   iter = 1

#   for j=1:ngroups
#     if j!=i 
#       pdensity[i,j] = v[iter]
#       iter+=1
#     end
#   end
# end

# pdensity = max.(pdensity,pdensity')

# A = sparse(stochastic_block_model(pdensity,repeat([round(Int,nnodes/ngroups)],ngroups)))
# A = largest_component(A)[1]

using Arpack

# # A,xy = spatial_network(50000, 2; degreedist=LogNormal(log(4),1.0),p=0);
# A,xy = spatial_network(50000, 2; degreedist=LogNormal(log(4),1.0),p=5);
# tdegs = vec(sum(A;dims=1))
# nnodes = lastindex(A,1)

# A,xy = adjust_graph(A,xy;factor=0.95,longedges=2.5,tdegs)
# nnz(A)/lastindex(A,1)


# p_vals = [5, 8, 10]
# factor_vals = [0.95, 0.99, 1.0, 1.05]
# longedges_vals = [0.75, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.25, 1.5]
# maxdeg = 50
# maxtrials = 200

# ps = vcat(Base.Iterators.product(p_vals,factor_vals,longedges_vals)...)



function get_stats(A::SparseMatrixCSC,xy::Matrix)
  nnodes = lastindex(A,1)
  #average degree 
  avgd = nnz(A)/lastindex(A,1)

  #spectral info 
  ls = eigs(A,nev=1)[1][1]

  #figure 
  f,xres,yres = ppr_sample_figure(A);
  plot!(f,yticks=[1e-3;1e-2;1e-1;1.0],ylims=(5e-4,1.1),
        xlims=(1,60000))

  #area 
  minx,miny = get_approximate_mins(xres,yres);
  extrapolate!(minx,miny,nnodes)
  area = approximate_auc(log10.(minx)./log10(nnodes/2), abs.(log10.(miny)))
  
  return vcat([avgd],[ls],[area]),f
end

dstpath = "/p/mnt/scratch/network-epi/study-testing/"

function single_param(p,factor,longedges,dstpath=dstpath)
  #base network
  A,xy = spatial_network(50000, 2; degreedist=LogNormal(log(4),1.0),p=5);
  tdegs = vec(sum(A;dims=1))
  
  data,fig = get_stats(A,xy)
  rounded = round.(data,digits=2)
  plot!(fig,title="iteration: 0\n avgd: $(rounded[1]), lambda1: $(rounded[2]), area: $(rounded[3])")
  Plots.savefig(fig,joinpath(dstpath,"$p-$factor-$longedges-0.png")) 

  #update file 
  open(joinpath(dstpath,"study-testing-$p-$factor-$longedges.txt"), "a") do io
    writedlm(io, [p factor longedges 0 data[1] data[2] data[3]])
  end

  #update graph 
  @showprogress for iter=1:maxtrials
    A,xy = adjust_graph(A,xy;factor=factor,longedges=longedges,tdegs)

    data,fig = get_stats(A,xy)
    rounded = round.(data,digits=2)
    plot!(fig,title="iteration: $iter\n avgd: $(rounded[1]), lambda1: $(rounded[2]), area: $(rounded[3])")

    Plots.savefig(fig,joinpath(dstpath,"$p-$factor-$longedges-$iter.png")) 
    #update file 
    open(joinpath(dstpath,"study-testing-$p-$factor-$longedges.txt"), "a") do io
      writedlm(io, [p factor longedges iter data[1] data[2] data[3]])
    end
  end
  return true 
end

#functions for measuring area
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


#parallelize this part 

# pmap(x->single_param(x[1],x[2],x[3]),ps)

# c = 0
# ind,m = 1,0.0
# for i=1:50
#   A,xy = adjust_graph(A,xy;factor=0.95,longedges=1.1,tdegs)
#   @show nnz(A)/lastindex(A,1)
#   c+=1

#   # f = plot(graphlines(A,xy),leg=false,alpha=0.1,markershape=:circle,
#   #   markerstrokewidth=0,markersize=2,framestyle=:none,margins=-5Measures.mm)
#   # A = loadGraph("study-11-2022-50.smat","input/graphs/")

#   f,xres,yres = ppr_sample_figure(A);
#   plot!(f,yticks=[1e-3;1e-2;1e-1;1.0],ylims=(5e-3,1.1),
#         xlims=(1,7500),title="iteration: $c")

#   #auc measure
#   minx,miny = get_approximate_mins(xres,yres);
#   extrapolate!(minx,miny,nnodes)
#   area = approximate_auc(log10.(minx)./log10(nnodes), abs.(log10.(miny))) 
#   @show area
#   if area>m
#     m = area
#     ind = i
#   end  
# end

# f

# l = eigs(A,nev=1)[1][1]

# function get_edge_lengths(A::SparseMatrixCSC,xy)
#   wv = zeros(Float64,nnz(A))
#   c = 1 
#   for i=1:size(A,1)
#     for nzi = nzrange(A,i)
#       j = rowvals(A)[nzi]
#       wv[c]= sqrt((xy[1,i]-xy[1,j])^2 + (xy[2,i]-xy[2,j])^2) / nnz(A)
#       c+=1
#     end
#   end  
#   return wv
# end

# vv = get_edge_lengths(A,xy)

# p = sortperm(vv,rev=true)
# vv[p]

# sum(vv.>0.75*mean(vv))
# nnz(A)


# #testing in parallel 
# #record average degree, lambda1, ppr area, and save fig 


# fnames = filter!(x->endswith(x,".txt"),readdir(dstpath))
# inds = [20;30:41;68:81]
# ind = 70

# fdata = readdlm(joinpath(dstpath,fnames[ind]))
# ind+=1
# scatter(fdata[:,5],fdata[:,7])

# fnames[30]

# fdata = readdlm(joinpath(dstpath,fnames[70]))

# # look at corresponding figures
# tmp = [1;6;11;21;51;101;201]
# fdata[tmp,7]


# fnames[70]



# gnames = getgnames("longrange-5-","input/graphs/")
# gname = gnames[end]
# A = loadGraph(gname,"input/graphs/")
# f,xres,yres = ppr_sample_figure(A)
# avgd = round(nnz(A)/lastindex(A,1),digits=2)
# title!(f,"PPR NCP: $(gname[1:end-5])\n Avgd: $(avgd)")





using Plots
using Measures

function spatial_graph_edges1(n::Integer,d::Integer;degreedist=LogNormal(log(4),1))
  xy = rand(d,n).- 0.5
  for col = 1:n
    xy[:,col]./=norm(xy[:,col])
  end
  T = BallTree(xy)
  # form the edges for sparse
  ei = Int[]
  ej = Int[]
  for i=1:n
    deg = min(ceil(Int,rand(degreedist)),n-1)
    idxs, dists = knn(T, xy[:,i], deg+1)
    for j in idxs
      if i != j
        push!(ei,i)
        push!(ej,j)
      end
    end
  end
  return xy, ei, ej
end
"""
n nodes in dimension d (random [0,1] box)
degreedist = degree distribution function
p = number of expected edges per node to add at ranodm (default = 0 )
"""
function spatial_network1(n::Integer, d::Integer; degreedist=LogNormal(log(3),1),
    p::Real=0.0)
  xy, ei, ej = spatial_graph_edges1(n, d;degreedist=degreedist)
  A = sparse(ei,ej,1,n,n)
  C = A
  if p > 0
    C .+= triu(sprand(Bool, n,n,p/n),1)
    fill(C.nzval,1)
  end
  return max.(C,C'), xy
end

A,xy = spatial_network1(50000,2,degreedist=LogNormal(log(5),1),p=0.5)
tdegs = vec(sum(A;dims=2))

# f = plot(graphlines(A,xy),leg=false,alpha=0.1,markershape=:circle,
#   markerstrokewidth=0,markersize=0,framestyle=:none,margins=-5Measures.mm)

A,xy = adjust_graph(A,xy,factor=0.95,longedges=3.0,tdegs=tdegs)

# f = plot(graphlines(A,xy),leg=false,alpha=0.1,markershape=:circle,
#   markerstrokewidth=0,markersize=2,framestyle=:none,margins=-5Measures.mm)

nnz(A)/lastindex(A,1)

f,_ = ppr_sample_figure(A,100)
f

l,v = eigs(A,nev=1)[1:2]

v = vec(v)
v./=sum(v)

p = sortperm(v,rev=true)
findfirst(cumsum(v[p]).>=0.999)

# A = loadGraph("study-11-2022-50.smat","input/graphs/")
kcore,_ = corenums(A)
kset = Set(sort(unique(kcore),rev=true)[1])

inds = [x for x=1:lastindex(A,1) if kcore[x] in kset]
inds = findall(vec(sum(A[inds,:];dims=1)).!=0)

sum(v[inds])

eigs(A,nev=5)[1]
entropy(v)/log(lastindex(A,1))



gnames = getgnames("study-11-2022-","input/graphs/")
push!(gnames,getgnames("noweak","input/graphs/")...)

using DelimitedFiles
A = loadGraph(gnames[1],"input/graphs/")
xy = readdlm(joinpath(mainDir,"../input/graphs/","$(gnames[end][1:end-5]).xy"))

get_stats


filter(x-> endswith(x,".xy"),readdir("/p/mnt/scratch/network-epi/input/graphs/"))


#################################################################
#original study graphs #remember to run in Julia1.5
using Random
using ProgressMeter

Random.seed!(10)
#noweak
A,xy = spatial_network(50000,2,degreedist=LogNormal(log(4),1),p=0.0)
writeSMAT("study-11-2023-noweak.smat", A;values=false)
writedlm("study-11-2023-noweak.xy",xy')

Random.seed!(10)
#other networks 
A,xy = spatial_network(50000,2,degreedist=LogNormal(log(4),1),p=0.5)
writeSMAT("study-11-2022-1.smat", A;values=false)
writedlm("study-11-2022-1.xy",xy')

@showprogress for iter=2:50
  A,xy = adjust_graph(A,xy;factor=0.95,longedges=1.0,tdegs)
  writeSMAT("study-11-2022-$iter.smat", A;values=false)
  writedlm("study-11-2022-$iter.xy",xy')  
end
###############################################################



gpath ="/p/mnt/scratch/network-epi/tmp-study/"
gnames = getgnames("noweak",gpath)
push!(gnames,getgnames("2022",gpath)...)

#doing degree histograms for david

figs = []
@showprogress for gname in gnames
  A = loadGraph(gname,gpath)
  # xy = readdlm(joinpath(gpath,"$(gname[1:end-5]).xy"))
  dvec = vec(sum(A;dims=2))
  f = histogram(dvec,
    bins=(0.9:1:250),leg=false,
    yscale=:log10,ylims = (0.8,7000),
    ylabel="frequency",xlabel="degree",
    title=gname)
  push!(figs,f)
end

figs[1]
figs[2]
figs[6]
figs[11]
figs[21]
figs[31]
figs[41]
figs[51]

#### spectral localization 

A = loadGraph(gnames[1],gpath)

ls,vs, _ = eigs(A,nev=10)
v = vs[:,1]
v./=sum(v)

p = sortperm(v,rev=true)
cs = cumsum(v[p])
findfirst(cs.>0.95)

plot(cs,leg=false,
  xscale=:log10,xlabel="number of nodes",
  ylabel="spectral localization")


#known sources of localization in networks 
#deepest kcore, star graph, and one-hop neighbors 
#kcore
kcore,_ = corenums(A)
kmax = maximum(kcore)

inds = findall(kcore.==kmax)
neighs = findall(vec(sum(A[inds,:],dims=1)).>0)
inds = union(inds,neighs)

sum(v[inds])

#star graph
dvec = vec(sum(A;dims=1))

inds = partialsortperm(dvec,1:10,rev=true)
neighs = findall(vec(sum(A[inds,:],dims=1)).>0)
inds = union(inds,neighs)

sum(v[inds])

localization_data = zeros(Float64,(51,2))
@showprogress for (iter,gname) in enumerate(gnames)
  A = loadGraph(gname,gpath)
  ls,vs, _ = eigs(A,nev=1)

  v = vs[:,1]
  v./=sum(v)

  p = sortperm(v,rev=true)
  cs = cumsum(v[p])

  xs,ys = get_approximate_mins(1:lastindex(cs),cs,5000)
  val = approximate_auc(xs,ys)/lastindex(cs)
  localization_data[iter,1] = iter
  localization_data[iter,2] = val
end
  
f = scatter(localization_data[:,1].-1,localization_data[:,2],label="Localization",
  color=1,markerstrokewidth=0,
  xlabel="Study Graph Iteration\n(0=noweak)",
  ylabel="Spectral Localization\n(AUC for λ₁,v₁)",
  title="study graphs spectral localization")


#spectral spacing 
neigs = 10
spacing_data = zeros(Float64,(51,neigs)) 
@showprogress for (iter,gname) in enumerate(gnames)
  A = loadGraph(gname,gpath)
  ls = first(eigs(A,nev=neigs))

  val = mean([ls[x]-ls[x+1] for x=1:lastindex(ls)-1])
  
  spacing_data[iter,1] = iter
  spacing_data[iter,2] = val
end


scatter(spacing_data[:,1].-1,spacing_data[:,2],label="Spacing",
  color=3,markerstrokewidth=0,
  ylabel="Spetral Spacing \nmean(Δ₁,Δ₂,..,Δ₉)",
  xlabel="Study Graph Iteration\n(0=noweak)",
  title="study graphs spectral spacing")

f = plot()
nmarkers = lastindex(localization_data,1)
marker_alphas = range(0.2,1,nmarkers)
for iter= 1:nmarkers
  scatter!([localization_data[iter,2]],[spacing_data[iter,2]],
  leg=false,color=2,
  markerstrokewidth=0,
  markershape=:circle,
  markeralpha=marker_alphas[iter])  
end
f

plot!(f,localization_data[:,2],spacing_data[:,2],
      linestyle=:dash,linewidth=0.25,
      linecolor=:black,
      xlabel="(λ₁,v₁) - Localization",
      ylabel="Average Spectral Spacing\n(Δ₁,…,Δ₉)",
      title="study graph evolution\n(light = 0, dark = 50)")



###testing for explicit geometry with a dense core
using Graphs
using SparseArrays

g,dists = Graphs.SimpleGraphs.euclidean_graph(xy,L=1.0,cutoff=0.2)

dists
A = sparse(g)

nnz(A)/lastindex(A,1)

dvec = vec(sum(A;dims=1))

using Plots
histogram(dvec)

#want to build a different type of graph with explicit geometry
#do euclidean_graph but need to be faster. 
#embed points in space 
#connect nodes of distance less than cutoff 
#nodes closer to "center" will be have higher degree


## starting over 
#= need functions that 
  - initializes a spatial network with a core bias
  - given a current set of points, reforms edges spatially
=#


#####initializing core bias. 
#maxdeg node in center 

degdist = LogNormal(log(4),1.0)
nnodes = 5000

#initialize pts 
xy = rand(2,nnodes)

#sample degs 
dvec = min.(ceil.(Int,rand(degdist,nnodes)),nnodes-1)
dmax,dind = findmax(dvec)

#### just sort distances by degrees 
p = sortperm([norm(xy[:,k].-0.5) for k=1:nnodes])
xy[:,p]
xy = xy[:,p]
sort!(dvec,rev=true)

#build spatial graph 

function spatial_network(xy::Matrix,dvec::Vector{Int})
  n = lastindex(xy,2)
  T = BallTree(xy)
  # form the edges for sparse
  ei = Int[]
  ej = Int[]
  for i=1:n
    deg = min(dvec[i],n-1)
    idxs, dists = knn(T, xy[:,i], deg+1)
    for j in idxs
      if i != j
        push!(ei,i)
        push!(ej,j)
      end
    end
  end
  A = sparse(ei,ej,1,n,n)
  A = max.(A,A')
  return A,xy
end

#=forces on nodes 
  - largest degree neighbor (if larger than our degree)
  - 
=#
#move towards biggest degree neighbor
#move towards core/periphery based on normalized degree

function adjust_positions(A,xy;factor=0.95,globalfactor=0.9)
  n = lastindex(A,1)
  center = [0.5;0.5]

  degs = vec(sum(A; dims=2))
  dmin,dmax = extrema(degs)

  centralization_factor = (degs.-dmin)./(dmax-dmin+0.25) 

  newxy = copy(xy)
  rowval = rowvals(A)
  for i=1:size(A,1)
    #moving towards center
    newxy[:,i] = centralization_factor[i]*newxy[:,i] + (1-centralization_factor[i]).*center
    bigj = i
    myd = degs[i]
    for nzi = nzrange(A,i)
      j = rowval[nzi]
      if degs[j] > degs[bigj]
        bigj = j
      end
    end
    if bigj != i # then we have a bigger neighbor
      newxy[:,i] = factor*newxy[:,i] + (1-factor)*newxy[:,bigj]  + 0.001*randn(2)
    end
  end
  return newxy
end

A,xy = spatial_network(xy,2*tdegs)
# xy = adjust_positions(A,xy,globalfactor=0.75)
# A,xy = spatial_network(xy,dvec)

scatter(xy[1,:],xy[2,:],
  markersize=1,markerstrokewidth=0,
  leg=false)
plot!(graphlines(A,xy),linewidth=0.05,leg=false)


nnz(A)/lastindex(A,1)

f,_ = ppr_sample_figure(A)
f

tdegs = vec(sum(A,dims=2))

A,xy = adjust_graph(A,xy;factor=0.95,longedges=0.95,tdegs)
f,_ = ppr_sample_figure(A)
f

nnz(A)/lastindex(A,1)


###### a (hopefully) final attempt at this
function spatial_graph_edges1(n::Integer,d::Integer;degreedist=LogNormal(log(4),1.1))
  xy = rand(d,n)
  dvec = min.(ceil.(Int,rand(degdist,nnodes)),nnodes-1)

  #### just sort distances by degrees 
  p = sortperm([norm(xy[:,k].-0.5) for k=1:nnodes])
  xy = xy[:,p]
  sort!(dvec,rev=true)

  T = BallTree(xy)
  # form the edges for sparse
  ei = Int[]
  ej = Int[]
  for i=1:n
    deg = min(ceil(Int,rand(degreedist)),n-1)
    idxs, dists = knn(T, xy[:,i], deg+1)
    for j in idxs
      if i != j
        push!(ei,i)
        push!(ej,j)
      end
    end
  end
  return xy, ei, ej
end
"""
n nodes in dimension d (random [0,1] box)
degreedist = degree distribution function
p = number of expected edges per node to add at ranodm (default = 0 )
"""
function spatial_network1(n::Integer, d::Integer; degreedist=LogNormal(log(3),1),
    p::Real=0.0)
  xy, ei, ej = spatial_graph_edges(n, d;degreedist=degreedist)
  A = sparse(ei,ej,1,n,n)
  C = A
  if p > 0
    C .+= triu(sprand(Bool, n,n,p/n),1)
    fill(C.nzval,1)
  end
  return max.(C,C'), xy
end


function plot_network(A,xy)
  f = scatter(xy[1,:],xy[2,:],
  markersize=1,markerstrokewidth=0,
  leg=false)
  plot!(f,graphlines(A,xy),linewidth=0.05,leg=false)
  return f 
end


A,xy = spatial_network1(3000,2,p=3.0,degreedist=LogNormal(log(4),1.0))
lastindex(largest_component(A)[1],1)
nnz(A)/lastindex(A,1)

# plot_network(A,xy)

f,_ = ppr_sample_figure(A)
f



xy = adjust_positions(A,xy)
plot_network(A,xy)

f,_ = ppr_sample_figure(A)
f

tdegs = vec(sum(A,dims=2))

A,xy = adjust_graph(A,xy;factor=0.95,longedges=2,tdegs)
f,_ = ppr_sample_figure(A)
f

nnz(A)/lastindex(A,1)

#for resetting when structure gets crazy 
A,xy = spatial_network(xy,tdegs)



pts = randn(2,500)
pts = mapslices(x->x/(2*norm(x))*rand(),pts,dims=1)

# scatter(pts[1,:],pts[2,:],leg=false,
#   markerstrokewidth=0,markersize=1)

T = BallTree(pts)

ei,ej = [],[]
@showprogress for node =1:lastindex(pts,2)
  idxs = inrange(T,pts[:,node],0.06);
  for k in idxs 
    push!(ei,node)
    push!(ej,k)
  end
end

A = sparse(ei,ej,1,lastindex(pts,2),lastindex(pts,2))
nnz(A)/lastindex(A,1)

plot_network(A,pts)
largest_component(A)[1]

histogram(vec(sum(A;dims=2)),leg=false)




###



A = loadGraph("study-24-100.smat","input/graphs/");
f,_ = epi_sample_figure(A,5e-1)

f