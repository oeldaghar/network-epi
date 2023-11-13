## ncp computations
using Distributed, Plots, Interact, DataFrames
using MatrixNetworks, SparseArrays, LinearAlgebra, Distributions, DelimitedFiles, ProgressMeter, Random, Plots, CSV


## ncp computation
include("diffusions.jl")
include("ncpplots1.jl")
# ENV["GKSwstype"] = "100"

# helper functions
function makeSparseMatrix(S::Vector{Set{Int}},nnodes::Int)
    ei,ej = Vector{Int}(),Vector{Int}()
    for (i,s) in enumerate(S)
        for k in s
            push!(ei,i)
            push!(ej,k)
        end
    end
    return sparse(ei,ej,ones(Int,length(ei)),length(S),nnodes)
end
function logscale_histogram(xs,ymin=1/length(xs);logx::Bool=false)
    h = float(StatsBase.fit(Histogram,xs,nbins=1+ceil(Int,maximum(abs.(xs)))))
    h.weights./=length(xs)

    #a little padding to make plots look better
    bins = [collect(h.edges[1]);maximum(h.edges[1])+1]
    w = [h.weights;0] .+ 0.5*ymin
    if logx
        f = plot(bins,w,seriestype=:barbins,label="",
                xlims=(1,2*maximum(xs)),xscale=:log10,
                ylims=(ymin,1),yscale=:log10)
    else
        f = plot(bins,w,seriestype=:barbins,label="",
                xlims=(0,maximum(xs)+2),ylims=(ymin,1),yscale=:log10)
    end
    return f
end #for plotting

function make_ncp(A::SparseMatrixCSC;
                    # giantcomponent::Bool=false,
                    get_sets::Bool = false,
                    rseed::Int=-1)

    @assert(isequal(size(A)...))

    A = SparseMatrixCSC{Int,Int}(A)

    if rseed<0
        rseed = Int(time_ns())
    end

    Random.seed!(rseed)
    nnodes = maximum(scomponents(A).sizes)
    if nnodes>4.5e4 #just a rule of thumb that ive noticed for number of sets that ncp finds
        a = 0.995
        mv = 5
    elseif 2e4<=nnodes<=4.5e4
        a = 0.99
        mv = 10
    elseif 1e4<=nnodes<2e4
        a = 0.95
        mv = 20
    else
        a = 0.8
        mv = 25
    end

    ncp,sets = DiffusionAlgorithms.serial_weighted_ncp(A;
        alpha=a, maxvisits=mv, get_sets=get_sets)
    return ncp,sets
end
# getgnames(gname,gpath) = sort(filter(x->occursin(lowercase(gname),lowercase(x)) && endswith(x,".smat"),readdir(gpath)),by=length)

function make_ncp(gname::String;gpath::String="",
            giantcomponent::Bool=true,
            get_sets::Bool = false,
            rseed::Int=-1,
            dst::String="")

    println("working on graph $gname")
    A = loadGraph(gname,gpath)
    if giantcomponent #do this by default, plots don't look right for disconnected graph
        A = largest_component(A)[1]
    end

    ncp,sets = make_ncp(A;get_sets=get_sets,rseed=rseed)
    #save ncp
    CSV.write(dst*"ncpinfo-$(gname[1:end-5]).txt",ncp)
    if get_sets#save the sets too
        writeSMAT(makeSparseMatrix(sets,size(A,1)),dst*"ncpmat-$gname")
    end
end

#plotting fcns 
function plot_ncp(gname::String;gpath::String="",dloc::String="",dst::String="")
    g = gname[1:end-5]
    A = loadGraph(gname,gpath)
    comps = scomponents(A)
    ncp = readdlm(dloc*"ncpinfo-$g.txt",',','\n')[2:end,:]
    #cols are seed,eps,size,cond,cut,volume_seed,volume_other,ComputeTime

    f = scatter(map(x->min(ncp[x,3],comps.sizes[comps.map[ncp[x,1]]]-ncp[x,3]),1:size(ncp,1)),ncp[1:end,4],
            xlabel = "size",ylabel = "conductance", xscale = :log10, yscale = :log10,
            xlims = (5,1.2*size(A,1)), ylims = (1e-5,1),
            title = "NCP - $g",label="")

    Plots.savefig(f,dst*"ncp-$g.png")
    return f
end

function plot_ncp(ncp::DataFrame,A::SparseMatrixCSC)

    f = scatter(map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1)),ncp.cond,scale=:log10,yscale=:log10,
        xlabel="size",ylabel="conductance",xlims=(5,1.2*size(A,1)),
        ylims=(1e-5,1),label="")

    f
end

function plot_hexbinncp(gname::String;gpath::String="",dloc::String="",dst::String="")
    g = gname[1:end-5]
    A = loadGraph(gname,gpath)
    ncp,headerinfo = readdlm(dloc*"ncpinfo-$(gname[1:end-5]).txt",',',header=true)
    ncp = DataFrame(ncp,vec(headerinfo))
    #cols are seed,eps,size,cond,cut,volume_seed,volume_other,ComputeTime

    f = myncpplot(map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1)),ncp.cond)
    Plots.plot!(xlabel="size",ylabel="conductance",xlims=(5,1.2*size(A,1)),
        ylims=(1e-5,1),title = "NCP - $g")

    Plots.savefig(f,dst*"ncphexbin-$g.png")
    return f
end

function plot_hexbinncp(ncp::DataFrame,A::SparseMatrixCSC)

    f = myncpplot(map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1)),ncp.cond)
    Plots.plot!(xlabel="size",ylabel="conductance",xlims=(5,1.2*size(A,1)),
        ylims=(1e-5,1),label="")

    f
end

function plot_hexbinncp_local(ncp::DataFrame,A::SparseMatrixCSC)
    """
        function for producing hexbin ncp plots locally (not on the server)
    """
    ncp.size = map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1))

    f = myncpplot_local(ncp.size,ncp.cond)
    Plots.plot!(xlabel="size",ylabel="conductance",xlims=(5,1.2*size(A,1)),
        ylims=(1e-4,1),label="")

    f
end

## fcns for coverage of ncp sets

function ncpmat_to_sets(ncpmat::SparseMatrixCSC)
    S = Vector{Set{Int}}(undef,size(ncpmat,1))
    @inbounds for i=1:size(ncpmat,1)
        S[i] = Set{Int}()
    end
    ei,ej,ev = findnz(ncpmat)
    @inbounds for k=1:length(ei)
        push!(S[ei[k]],ej[k])
    end
    S
end

function ncpcoverage(gname::String;dloc::String="",dst::String="")
    #load ncp and sets
    ncp,headerinfo = readdlm(dloc*"ncpinfo-$(gname[1:end-5]).txt",',',header=true)
    ncp = DataFrame(ncp,vec(headerinfo))
    ncpmat = MatrixNetworks.readSMAT(dloc*"ncpmat-$gname")
    sets = ncpmat_to_sets(ncpmat)

    nnodes = size(ncpmat,2)
    #modifying to always use the smaller sets
    for i=1:length(sets)
        s = sets[i]
        if length(s)>nnodes/2
            sets[i] = Set(deleteat!(collect(1:nnodes),sort(collect(s))))
        end
    end

    p = sortperm(ncp.cond)
    cvals = sort(unique(ncp.cond))
    S = Set{Int}()
    j = 1
    v = cvals[j]
    vals = zeros(Float64,length(unique(ncp.cond)))
    for (i,s) in enumerate(sets[p])
        if ncp.cond[p[i]]>v
            vals[j] = length(S)
            j+=1
            v = cvals[j]
        end
        union!(S,s)
    end

    ymin = floor(log10(cvals[1]))-1
    f = scatter(vals[1:end-1]./nnodes,cvals[1:end-1],
        ylims=(1e-4,1.1),xlims=(0,1.05),yscale=:log10,label="",
        ylabel="conductance",xlabel="coverage",
        title="$(gname[1:end-5])\ncoverage vs conductance")

    Plots.savefig(f,dst*"ncpcoverage-$(gname[1:end-5]).png")
    return f
end

function ncpcoverage(A::SparseMatrixCSC;dst::String="")
    #get ncp and sets
    ncp,sets = make_ncp(A,giantcomponent=true,get_sets=true)

    nnodes = size(A,2)
    #modifying to always use the smaller sets
    for i=1:length(sets)
        s = sets[i]
        if length(s)>nnodes/2
            sets[i] = Set(deleteat!(collect(1:nnodes),sort(collect(s))))
        end
    end

    p = sortperm(ncp.cond)
    cvals = sort(unique(ncp.cond))
    S = Set{Int}()
    j = 1
    v = cvals[j]
    vals = zeros(Float64,length(unique(ncp.cond)))
    for (i,s) in enumerate(sets[p])
        if ncp.cond[p[i]]>v
            vals[j] = length(S)
            j+=1
            v = cvals[j]
        end
        union!(S,s)
    end

    ymin = floor(log10(cvals[1]))-1
    f = scatter(vals[1:end-1]./nnodes,cvals[1:end-1],
        ylims=(1e-4,1.1),xlims=(0,1.05),yscale=:log10,label="",
        ylabel="conductance",xlabel="coverage")
        # ,title="$(gname[1:end-5])\ncoverage vs conductance")

    # Plots.savefig(f,dst*"ncpcoverage-$(gname[1:end-5]).png")
    return f
end

function ncpcoverage(ncp::DataFrame,sets::Vector{Set{Int}},A::SparseMatrixCSC)
    nnodes = size(A,2)
    #modifying to always use the smaller sets
    for i=1:length(sets)
        s = sets[i]
        if length(s)>nnodes/2
            sets[i] = Set(deleteat!(collect(1:nnodes),sort(collect(s))))
        end
    end

    p = sortperm(ncp.cond)
    cvals = sort(unique(ncp.cond))
    S = Set{Int}()
    j = 1
    v = cvals[j]
    vals = zeros(Float64,length(unique(ncp.cond)))
    for (i,s) in enumerate(sets[p])
        if ncp.cond[p[i]]>v
            vals[j] = length(S)
            j+=1
            v = cvals[j]
        end
        union!(S,s)
    end

    ymin = floor(log10(cvals[1]))-1
    f = scatter(vals[1:end-1]./nnodes,cvals[1:end-1],
        ylims=(1e-4,1.1),xlims=(0,1.05),yscale=:log10,label="",
        ylabel="conductance",xlabel="coverage")
    return f
end

function make_ncp_plots(gname::String;gpath::String="",dst::String="")
    A = loadGraph(gname,gpath)
    # A = largest_component(A)[1]
    davg = nnz(A)/size(A,1)
    ncp,sets = make_ncp(A,get_sets=true)
    f = plot_ncp(ncp,A)
    title!(f,"NCP - $(gname[1:end-5])\ndavg: $(round(davg,digits=2)), nnodes: $(size(A,1))")
    Plots.savefig(dst*"ncp-$(gname[1:end-5]).png")

    f = ncpcoverage(ncp,sets,A)
    title!(f,"Coverage - $(gname[1:end-5])\ndavg: $(round(davg,digits=2)), nnodes: $(size(A,1))")
    Plots.savefig(dst*"ncpcoverage-$(gname[1:end-5]).png")

    #save raw data
    CSV.write(dst*"ncpinfo-$(gname[1:end-5]).txt",ncp)
    writeSMAT(makeSparseMatrix(sets,size(A,1)),dst*"ncpmat-$gname")
end

function make_ncp_plots(ncp::DataFrame,sets::Vector{Set{Int}},A::SparseMatrixCSC;gname::String="",dst::String="")
    A = largest_component(A)[1]
    A = max.(A,A')
    davg = nnz(A)/size(A,1)
    f = plot_ncp(ncp,A)
    title!(f,"NCP - $(gname[1:end-5])\ndavg: $(round(davg,digits=2)), nnodes: $(size(A,1))")
    Plots.savefig(dst*"ncp-$(gname[1:end-5]).png")

    f = ncpcoverage(ncp,sets,A)
    title!(f,"Coverage - $(gname[1:end-5])\ndavg: $(round(davg,digits=2)), nnodes: $(size(A,1))")
    Plots.savefig(dst*"ncpcoverage-$(gname[1:end-5]).png")
end

function make_hexbinncp_plots(gname::String;gpath::String="",dst::String="")
    A = loadGraph(gname,gpath)
    A = largest_component(A)[1]

    davg = nnz(A)/size(A,1)
    ncp,sets = make_ncp(A,get_sets=true)

    #coverage
    f = ncpcoverage(ncp,sets,A)
    title!(f,"Coverage - $(gname[1:end-5])\ndavg: $(round(davg,digits=2)), nnodes: $(size(A,1))")
    Plots.savefig(f,dst*"ncpcoverage-$(gname[1:end-5]).png")

    #ncp
    f = plot_hexbinncp(ncp,A)
    title!(f,"NCP - $(gname[1:end-5])\ndavg: $(round(davg,digits=2)), nnodes: $(size(A,1))")
    Plots.savefig(f,dst*"ncphexbin-$(gname[1:end-5]).png")

    #save raw data
    CSV.write(dst*"ncpinfo-$(gname[1:end-5]).txt",ncp)
    writeSMAT(makeSparseMatrix(sets,size(A,1)),dst*"ncpmat-$gname")
    return f
end

function make_hexbinncp_plots(ncp::DataFrame,sets::Vector{Set{Int}},A::SparseMatrixCSC;gname::String="",dst::String="")
    A = largest_component(A)[1]
    A = max.(A,A')
    davg = nnz(A)/size(A,1)

    # f = ncpcoverage(ncp,sets,A)
    # title!(f,"Coverage - $(gname[1:end-5])\ndavg: $(round(davg,digits=2)), nnodes: $(size(A,1))")
    # Plots.savefig(dst*"ncpcoverage-$(gname[1:end-5]).png")

    f = plot_hexbinncp(ncp,A)
    title!(f,"NCP - $(gname[1:end-5])\ndavg: $(round(davg,digits=2)), nnodes: $(size(A,1))")
    # Plots.savefig(dst*"ncphexbin-$(gname[1:end-5]).png")
    return f
end

function load_ncpdata(gname::String,dloc::String="")
    if endswith(gname,".smat")
        g = gname[1:end-5]
    else
        g = gname
    end

    ncp,headerinfo = readdlm(dloc*"ncpinfo-$g.txt",',',header=true)
    ncp = DataFrame(ncp,vec(headerinfo))
    ncpmat = MatrixNetworks.readSMAT(dloc*"ncpmat-$g.smat")
    sets = ncpmat_to_sets(ncpmat)
    return ncp,ncpmat,sets
end

##for local hexbin plotting.
#fill_z appears to have undergone a change after v1.5.1 
#this is copied and modified from
#https://raw.githubusercontent.com/RalphAS/HexBinPlots.jl/master/src/HexBinPlots.jl

using Hexagons
function xy2inds_(x::AbstractArray,y::AbstractArray,
    bins::Union{NTuple{1,Int},NTuple{2,Int},Int})

    if length(bins) == 2
        xbins, ybins = bins
    else
        xbins, ybins = (bins...,bins...)
    end
    xmin, xmax = extrema(x)
    ymin, ymax = extrema(y)
    xspan, yspan = xmax - xmin, ymax - ymin
    xsize, ysize = xspan / xbins, yspan / ybins
    x0, y0 = xmin - xspan / 2,ymin - yspan / 2

    inds = Dict{(Tuple{Int, Int}), Vector{Int}}()
    cnts = Dict{(Tuple{Int, Int}), Int}()

    for i in eachindex(x)
        h = convert(HexagonOffsetOddR,
                    cube_round(x[i] - x0 + xsize, y[i] - y0 + ysize, xsize, ysize))
        idx = (h.q, h.r)

        cnts[idx] = 1 + get(cnts,idx,0)
        inds[idx] = push!(get(inds,idx,Vector{Int}()),i)
    end
    inds,xsize,ysize,x0,y0
end
function inds2xy_(inds::Dict{S,T},xsize, ysize, x0, y0) where {S,T}
    nhex = length(inds)
    xh = zeros(nhex)
    yh = zeros(nhex)
    vh = zeros(Int,nhex)
    k = 0
    for (idx, s) in inds
        k += 1
        xx,yy = Hexagons.center(HexagonOffsetOddR(idx[1], idx[2]),
        xsize, ysize, x0, y0)
        xh[k] = xx
        yh[k] = yy
        vh[k] = length(s)
    end
    xh,yh,vh
end

#the following functions use the ACL NCP and weigh the bins using the sets 
#from the ncp computation and the snodes specified from the diffusion data 

#weighs each bin by weighing each fraction of susceptible nodes in all sets in bin
# so binweight = mean_{nodes in bin} (number of times node was susceptible in all trials)
function plot_missed_sets(ncp::DataFrame,ncpmat::SparseMatrixCSC,snodes::Vector,sets::Vector{Set{Int}};
            nbins::Int=100,ntrials::Int=1)

    @assert(minimum(ncp.cond)>0,"graph not connected. this will throw an error in the hexbin plotting")
    @assert(maximum(snodes)<=ntrials,"ntrials < max # of times a node is susceptible")
    #larger values means we missed more nodes in that set
    missedsets = ncpmat*snodes./sum(ncpmat;dims=2)
    #for aggregating over different diffusions
    missedsets./=ntrials

    x,y = log10.(ncp.size),log10.(ncp.cond)
    inds,xsize,ysize,x0,y0 = xy2inds_(x,y,nbins)
    xh,yh,vh = inds2xy_(inds,xsize,ysize,x0,y0)
    hexhist = HexBinPlots.HexHistogram(xh,yh,vh,xsize,ysize,false)
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

    #normalizing bins
    wh = Vector{Float64}()
    # sets = ncpmat_to_sets(ncpmat)
    for (idx,s) in inds
        #fraction of susceptible nodes in bin
        nodes = collect(union(sets[s]...))
        push!(wh,sum(snodes[nodes])/length(nodes))

        #mean_{sets in bin}(frac of set susceptible)
        # push!(wh,sum(missedsets[s])/length(s))
    end

    #couldn't get fill_z to work. might have undergone a change since v1.5.1
    #so will jerry-rig this to plot one by one.
    f = Plots.plot(10.0.^xsh[1:8],10.0.^ysh[1:8],fill_z=wh[1],linecolor=nothing,
        seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=true,
        ylims=(1e-2,1),xlims=(3,2*maximum(x)),leg=false,
        clims=(0,1))
    for k=2:length(wh)
        winds = (k-1)*8+1:8*k
        Plots.plot!(10.0.^xsh[winds],10.0.^ysh[winds],fill_z=wh[k],linecolor=nothing,
            seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=true,
            ylims=(1e-2,1),xlims=(3,2*maximum(x)),leg=false,
            clims=(0,1))
    end
    ylims!(f,(1e-4,1))
    return f
end

function plot_missed_sets(gname::String;ncploc::String="",sloc::String="",
            nbins::Int=100,ntrials::Int=1)
    
    if endswith(gname,".smat")
        g = gname[1:end-5]
    else
        g = gname
    end

    ncp,headerinfo = readdlm(ncploc*"ncpinfo-$g.txt",',',header=true)
    ncp = DataFrame(ncp,vec(headerinfo))

    ncpmat = MatrixNetworks.readSMAT(ncploc*"ncpmat-$g.smat")
    sets = ncpmat_to_sets(ncpmat)

    #load snodes data
    snodes = readdlm(sloc*"scounts-$g-seir-0.05-0.05.txt",',',Int)
    if size(snodes,2)>1
        snodes = snodes[:,1] #using the first column #no quarantining
    end
    
    return plot_missed_sets(ncp,ncpmat,snodes,sets;nbins=nbins,ntrials=ntrials)
end


#weighs each bin by weighing each set by fraction of set susceptible then averaging over sets in the bin
# so binweight = sum_{sets in bin} (frac of set that is infected)
#old function for plotting in v1.5.1
# function plot_missed_sets1(ncp::DataFrame,ncpmat::SparseMatrixCSC,snodes,sets;
#             nbins::Int=100,ntrials::Int=1,color=:inferno)

#     @assert(maximum(snodes)<=ntrials,"ntrials < max # of times a node is susceptible")
#     #larger values means we missed more nodes in that set
#     missedsets = ncpmat*snodes./sum(ncpmat;dims=2)
#     #for aggregating over different diffusions
#     missedsets./=ntrials

#     x,y = log10.(ncp.size),log10.(ncp.cond)
#     inds,xsize,ysize,x0,y0 = xy2inds_(x,y,nbins)
#     xh,yh,vh = inds2xy_(inds,xsize,ysize,x0,y0)
#     hexhist = HexBinPlots.HexHistogram(xh,yh,vh,xsize,ysize,false)
#     h,vh = HexBinPlots.make_shapes(hexhist)

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

#     #normalizing bins
#     wh = Vector{Float64}()
#     for (idx,s) in inds
#         #fraction of susceptible nodes in bin
#         nodes = collect(union(sets[s]...))
#         push!(wh,(sum(snodes[nodes])/length(nodes))/ntrials)

#         #mean_{sets in bin} (frac of set susceptible)
#         # push!(wh,sum(missedsets[s])/length(s))
#     end

#     #couldn't get fill_z to work. might have undergone a change since v1.5.1
#     #so will jerry-rig this to plot one by one.
#     f = Plots.plot(10.0.^xsh[1:8],10.0.^ysh[1:8],fill_z=wh[1],linecolor=nothing,
#         seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=true,
#         ylims=(1e-2,1),xlims=(3,2*maximum(x)),leg=false,linewidth=0,
#         clims=(0,1),color=color)
#     for k=2:length(wh)
#         winds = (k-1)*8+1:8*k
#         Plots.plot!(f,10.0.^xsh[winds],10.0.^ysh[winds],fill_z=wh[k],linecolor=nothing,
#             seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=true,
#             ylims=(1e-2,1),xlims=(3,2*maximum(ncp.size)),leg=false,linewidth=0,
#             clims=(0,1),color=color)
#     end

#     #making colors looks better 
#     for k=1:length(wh)
#         winds = (k-1)*8+1:8*k
#         Plots.plot!(f,10.0.^xsh[winds],10.0.^ysh[winds],fill_z=wh[k],linecolor=nothing,
#             seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=true,
#             ylims=(1e-2,1),xlims=(3,2*maximum(ncp.size)),leg=false,linewidth=0,
#             clims=(0,1),color=color)
#     end

#     ylims!(f,(1e-4,1))
#     return f
# end

# updated function for v1.8.0 plotting
function plot_missed_sets1!(f,ncp::DataFrame,ncpmat::SparseMatrixCSC,snodes,sets;
            nbins::Int=100,ntrials::Int=1,color=:inferno)

    @assert(maximum(snodes)<=ntrials,"ntrials < max # of times a node is susceptible")
    #larger values means we missed more nodes in that set
    missedsets = ncpmat*snodes./sum(ncpmat;dims=2)
    #for aggregating over different diffusions
    missedsets./=ntrials

    x,y = log10.(ncp.size),log10.(ncp.cond)
    inds,xsize,ysize,x0,y0 = xy2inds_(x,y,nbins)
    xh,yh,vh = inds2xy_(inds,xsize,ysize,x0,y0)
    hexhist = HexBinPlots.HexHistogram(xh,yh,vh,xsize,ysize,false)
    h,vh = HexBinPlots.make_shapes(hexhist)


    #new
    xsh = Vector{Float64}()
    ysh = Vector{Float64}()
    vsh = Vector{Float64}()
    for k in eachindex(h)
      append!(xsh,h[k].x)
      push!(xsh,h[k].x[1])
      push!(xsh,NaN)
      append!(ysh,h[k].y)
      push!(ysh,h[k].y[1])
      push!(ysh,NaN)
      for i in 1:length(h[k].x)+2
        push!(vsh, vh[k])
      end 
    end
    
    #normalizing bins
    wh = Vector{Float64}()
    for (idx,s) in inds
        #fraction of susceptible nodes in bin
        nodes = collect(union(sets[s]...))
        bin_val = (sum(snodes[nodes])/length(nodes))/ntrials 
        for k=1:8
            push!(wh,bin_val)
        end
    end

    Plots.plot!(f,10.0.^xsh,10.0.^ysh,fill_z=wh,linecolor=nothing,
        seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=false,color=color)
        # xlims=xlims,ylims=ylims,color=color)
    
    Plots.plot!(f,10.0.^xsh,10.0.^ysh,fill_z=wh,linecolor=nothing,
        seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=false,color=color)
        # xlims=xlims,ylims=ylims)
    
    ylims!(f,(1e-4,1))
    return f
end


function plot_missed_sets1(ncp::DataFrame,ncpmat::SparseMatrixCSC,snodes,sets;
            nbins::Int=100,ntrials::Int=1,color=:inferno)

    @assert(maximum(snodes)<=ntrials,"ntrials < max # of times a node is susceptible")
    #larger values means we missed more nodes in that set
    missedsets = ncpmat*snodes./sum(ncpmat;dims=2)
    #for aggregating over different diffusions
    missedsets./=ntrials

    x,y = log10.(ncp.size),log10.(ncp.cond)
    inds,xsize,ysize,x0,y0 = xy2inds_(x,y,nbins)
    xh,yh,vh = inds2xy_(inds,xsize,ysize,x0,y0)
    hexhist = HexBinPlots.HexHistogram(xh,yh,vh,xsize,ysize,false)
    h,vh = HexBinPlots.make_shapes(hexhist)


    #new
    xsh = Vector{Float64}()
    ysh = Vector{Float64}()
    vsh = Vector{Float64}()
    for k in eachindex(h)
      append!(xsh,h[k].x)
      push!(xsh,h[k].x[1])
      push!(xsh,NaN)
      append!(ysh,h[k].y)
      push!(ysh,h[k].y[1])
      push!(ysh,NaN)
      for i in 1:length(h[k].x)+2
        push!(vsh, vh[k])
      end 
    end
    
    #normalizing bins
    wh = Vector{Float64}()
    for (idx,s) in inds
        #fraction of susceptible nodes in bin
        nodes = collect(union(sets[s]...))
        bin_val = (sum(snodes[nodes])/length(nodes))/ntrials 
        for k=1:8
            push!(wh,bin_val)
        end
    end

    f = Plots.plot(10.0.^xsh,10.0.^ysh,fill_z=wh,linecolor=nothing,
        seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=false,color=color)
        # xlims=xlims,ylims=ylims,color=color)
    
    Plots.plot!(f,10.0.^xsh,10.0.^ysh,fill_z=wh,linecolor=nothing,
        seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=false,color=color)
        # xlims=xlims,ylims=ylims)
    
    ylims!(f,(1e-4,1))
    return f
    
end

function plot_missed_sets1(gname::String;ncploc::String="",sloc::String="",
            nbins::Int=100,ntrials::Int=1,model::String="seir",
            betas::Vector{Float64}=vcat(1e-3:1e-3:1e-2, 2e-2:1e-2:1e-1),
            qpercents::Vector{Int}=0:15,
            color=:inferno)
    
    if endswith(gname,".smat")
        g = gname[1:end-5]
    else
        g = gname
    end

    ncp,headerinfo = readdlm(ncploc*"ncpinfo-$g.txt",',',header=true)
    ncp = DataFrame(ncp,vec(headerinfo))

    ncpmat = MatrixNetworks.readSMAT(ncploc*"ncpmat-$g.smat")
    sets = ncpmat_to_sets(ncpmat)
    figs = Dict{Tuple,Plots.Plot}()
    @showprogress for beta in betas
        #load snodes data
        snodes = readdlm(sloc*"scounts-$g-$model-$beta-0.05.txt",',',Int)
        for qpercent in qpercents #1:size(snodes)[2]
            snodes_col = snodes[:,qpercent+1] 
            avg_inf_frac = round(1-( sum(snodes_col) / (ntrials*lastindex(snodes_col)) ),digits=2)*100
            f =  plot_missed_sets1(ncp,ncpmat,snodes_col,sets;nbins=nbins,ntrials=ntrials,color=color)
            Plots.plot!(f,title="$gname - $model($beta,0.05,$(qpercent) Qpercent)\nAvg. Infected Fraction: $avg_inf_frac%",colorbar_title="Fraction Nodes Susceptible")
            figs[(beta,qpercent)] = f #(beta,qpercent): missed sets fig
        end
    end
    return figs 
end

function plot_missed_sets1(gname::String;ncploc::String="",sloc::String="",
            nbins::Int=100,ntrials::Int=1,model::String="seir",
            betas::Vector{Float64}=vcat(1e-3:1e-3:1e-2, 2e-2:1e-2:1e-1),
            qpercent::Int=0,color=:inferno)

    scol = qpercent+1
    
    if endswith(gname,".smat")
        g = gname[1:end-5]
    else
        g = gname
    end

    ncp,headerinfo = readdlm(ncploc*"ncpinfo-$g.txt",',',header=true)
    ncp = DataFrame(ncp,vec(headerinfo))

    ncpmat = MatrixNetworks.readSMAT(ncploc*"ncpmat-$g.smat")
    sets = ncpmat_to_sets(ncpmat)
    figs = Dict{Tuple,Plots.Plot}()
    @showprogress for beta in betas
        #load snodes data
        snodes = readdlm(sloc*"scounts-$g-$model-$beta-0.05.txt",',',Int)
        snodes_col = snodes[:,scol] 

        avg_inf_frac = round(1-( sum(snodes_col) / (ntrials*lastindex(snodes_col)) ),digits=2)*100
        f =  plot_missed_sets1(ncp,ncpmat,snodes_col,sets;nbins=nbins,ntrials=ntrials,color=color)
        Plots.plot!(f,title="$gname - $model($beta,0.05,$(qpercent) Qpercent)\nAvg. Infected Fraction: $avg_inf_frac%",colorbar_title="Fraction Nodes Susceptible")
        figs[(beta,qpercent)] = f #(beta,qpercent): missed sets fig
    end
    return figs 
end

function plot_missed_sets1!(f,gname::String;ncploc::String="",sloc::String="",
            nbins::Int=100,ntrials::Int=1,model::String="seir",
            betas::Vector{Float64}=vcat(1e-3:1e-3:1e-2, 2e-2:1e-2:1e-1),
            qpercent::Int=0,color=:inferno)

    scol = qpercent+1
    
    if endswith(gname,".smat")
        g = gname[1:end-5]
    else
        g = gname
    end

    ncp,headerinfo = readdlm(ncploc*"ncpinfo-$g.txt",',',header=true)
    ncp = DataFrame(ncp,vec(headerinfo))

    ncpmat = MatrixNetworks.readSMAT(ncploc*"ncpmat-$g.smat")
    sets = ncpmat_to_sets(ncpmat)
    figs = Dict{Tuple,Plots.Plot}()

    
    @showprogress for beta in betas
        f_loc = deepcopy(f)
        #load snodes data
        snodes = readdlm(sloc*"scounts-$g-$model-$beta-0.05.txt",',',Int)
        snodes_col = snodes[:,scol] 

        avg_inf_frac = round(1-( sum(snodes_col) / (ntrials*lastindex(snodes_col)) ),digits=2)*100
        f_loc = plot_missed_sets1!(f_loc,ncp,ncpmat,snodes_col,sets;nbins=nbins,ntrials=ntrials,color=color)
        Plots.plot!(f_loc,title="$gname - $model($beta,0.05,$(qpercent) Qpercent)\nAvg. Infected Fraction: $avg_inf_frac%",colorbar_title="Fraction Nodes Susceptible")
        figs[(beta,qpercent)] = f_loc #(beta,qpercent): missed sets fig
    end
    return figs 
end


#weighs each bin by weighing each set by fraction of set susceptible then averaging over sets in the bin
# so binweight = sum_{sets in bin} (frac of set that is infected)
function normalize_bins(h,vsh)
    centroid_to_inds = Dict()
    for k in eachindex(h)
        centroid = mean(h[k].x) #x-component only
        if centroid in keys(centroid_to_inds)
            append!(centroid_to_inds[centroid],collect(8*(k-1)+1:8*(k-1)+8))
        else
            centroid_to_inds[centroid] = Vector{Int}(collect(8*(k-1)+1:8*(k-1)+8))
        end
    end

    new_vals = zeros(Float64,lastindex(vsh))
    for k in keys(centroid_to_inds)
        inds = centroid_to_inds[k]
        new_vals[inds] .= vsh[inds]/sum(vsh[inds]) #normalize
    end
    return new_vals
end

function plot_missed_sets2(ncp::DataFrame,ncpmat::SparseMatrixCSC,snodes;
            nbins::Int=100,ntrials::Int=1,color=:inferno)

    @assert(maximum(snodes)<=ntrials,"ntrials < max # of times a node is susceptible")
    #larger values means we missed more nodes in that set
    missedsets = ncpmat*snodes./sum(ncpmat;dims=2)
    #for aggregating over different diffusions
    missedsets./=ntrials

    x,y = log10.(ncp.size),log10.(ncp.cond)
    inds,xsize,ysize,x0,y0 = xy2inds_(x,y,nbins)
    xh,yh,vh = inds2xy_(inds,xsize,ysize,x0,y0)
    hexhist = HexBinPlots.HexHistogram(xh,yh,vh,xsize,ysize,false)
    h,vh = HexBinPlots.make_shapes(hexhist)

    xsh = Vector{Float64}()
    ysh = Vector{Float64}()
    vsh = Vector{Float64}()
    for k in eachindex(h)
        append!(xsh,h[k].x)
        push!(xsh,h[k].x[1])
        push!(xsh,NaN)
        append!(ysh,h[k].y)
        push!(ysh,h[k].y[1])
        push!(ysh,NaN)
        for i in 1:length(h[k].x)+2
            push!(vsh, vh[k])
        end 
    end

    #normalizing bins
    wh = Vector{Float64}()
    # sets = ncpmat_to_sets(ncpmat)
    for (idx,s) in inds
        #fraction of susceptible nodes in bin
        # nodes = collect(union(sets[s]...))
        # push!(wh,sum(snodes[nodes])/length(nodes))

        #mean_{sets in bin} (frac of set susceptible)
        for i in 1:8
            push!(wh,sum(missedsets[s])/length(s))
        end
    end

    #normalizing bins by maximum 
    bin_weights = normalize_bins(h,wh)
    

    f = Plots.plot(10.0.^xsh,10.0.^ysh,fill_z=log10.(bin_weights.+1)./log10(2),
            seriestype=:shape,
            xscale=:log10,yscale=:log10,
            label="",
            xlims=(3,10.0^maximum(x)),ylims=(min(1e-4,10.0^minimum(y)/5),1.1),
            linecolor=nothing, colorbar=true, c=color, clims=(0,1)
        )
    return f
end

function plot_missed_sets2(gname::String;ncploc::String="",sloc::String="",
            nbins::Int=100,ntrials::Int=1,model::String="seir",
            betas::Vector{Float64}=vcat(1e-3:1e-3:1e-2, 2e-2:1e-2:1e-1))
    
    if endswith(gname,".smat")
        g = gname[1:end-5]
    else
        g = gname
    end

    ncp,headerinfo = readdlm(ncploc*"ncpinfo-$g.txt",',',header=true)
    ncp = DataFrame(ncp,vec(headerinfo))

    ncpmat = MatrixNetworks.readSMAT(ncploc*"ncpmat-$g.smat")
    figs = Dict()
    for beta in betas
        #load snodes data
        snodes = readdlm(sloc*"scounts-$g-$model-$beta-0.05.txt",',',Int)
        if size(snodes,2)>1
            snodes = snodes[:,1] #using the first column #no quarantining
        end
        figs[beta] = plot_missed_sets2(ncp,ncpmat,snodes;nbins=nbins,ntrials=ntrials)
    end
    return figs 
end
#=
function delta_hexbin(x1,y1,x2,y2;nbins=100)
    """given ncp info, (x1,y1) and (x2,y2), display ncp1 and info on ncp2 but not ncp1"""
  
    #merge info, bin data, then color bins based on ncp1
  
    x_merge,y_merge = vcat(x1,x2),vcat(y1,y2)
    xbounds = extrema(x_merge)
    ybounds = extrema(y_merge)  
    
    #bin data 
    hexhist = fit(HexBinPlots.HexHistogram,log10.(x_merge),log10.(y_merge),nbins,
      xlims = log10.(xbounds),ylims=log10.(ybounds))
  
    #sort ncp1 data for easy lookup
  
    h,vh = HexBinPlots.make_shapes(hexhist)
    xsh = Vector{Float64}()
    ysh = Vector{Float64}()
    vsh = Vector{Float64}()
    for k in eachindex(h)
      append!(xsh,h[k].x)
      push!(xsh,h[k].x[1])
      push!(xsh,NaN)
      append!(ysh,h[k].y)
      push!(ysh,h[k].y[1])
      push!(ysh,NaN)
      for i in 1:length(h[k].x)+2 
        #check if ncp1 contributed to bin 
        push!(vsh, vh[k])
      end 
    end
  
    Plots.plot!(f,10.0.^xsh,10.0.^ysh,fill_z=log10.(vsh.+1),linecolor=nothing,
          seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=false,
          xlims=xlims,ylims=ylims)
    f
  
  end
  =#