## ncp computations
using Distributed, Plots, Interact, DataFrames
using MatrixNetworks, SparseArrays, LinearAlgebra, Distributions, DelimitedFiles, ProgressMeter, Random, Plots, CSV

#allows for execution from command line as well as an ide so files can be run modularlly
mainDir = joinpath(split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))])

## ncp computation
include(joinpath(mainDir,"code","ncp","diffusions.jl"))
include(joinpath(mainDir,"code","ncp","ncpplots1.jl"))
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
