#for converting flow data into a contact network

data_dir = "/homes/oeldagha/COVID19USFlows/"
data_dir = "/homes/oeldagha/COVID19USFlows/ct2ct/"
dst = "/p/mnt/scratch/network-epi/input/graphs/"
fnames = filter(c->endswith(c,".csv"),readdir(data_dir))

year_fnames = filter(c->occursin("_2020_",c),fnames)

using DataFrames
using DelimitedFiles
using SparseArrays


function load_flowdata(fnames::Vector{String})
    #return variables
    ei,ej,ev = Vector{Int64}(),Vector{Int64}(),Vector{Int64}()
    nodes = Set{Tuple{Int,Float64,Float64}}() #geoid, longitude, latitude

    for fname in fnames 
        println("working on $fname")
        f = open(fname)
        #process header 
        header = readline(f)
        while !eof(f)
            line = readline(f)
            parts = split(line,',')
            parts = split(line,',')
            
            geoid_o = parse(Int64, parts[1])
            geoid_d = parse(Int64, parts[2])
            lng_o = parse(Float64, parts[3])
            lat_o = parse(Float64, parts[4])
            lng_d = parse(Float64, parts[5])
            lat_d = parse(Float64, parts[6])
            visitor_counts = parse(Int64, parts[8])

            push!(ei,geoid_o)
            push!(ej,geoid_d)
            
            #keep weights for testing
            push!(ev,visitor_counts)

            push!(nodes,(geoid_o,lng_o,lat_o))
            push!(nodes,(geoid_d,lng_d,lat_d))
        end
    end

    #geoid -> new node id 
    nodemap = Dict{Int,Int}(map(x->(x[2][1],x[1]),enumerate(nodes)))
    nnodes = length(nodemap)

    ei = map(x->nodemap[x],ei)
    ej = map(x->nodemap[x],ej)
    return ei,ej,ev,nodes
    # A = sparse(ei,ej,ev,n,n)
    # return A,nodes,nodemap 
end
year_fnames
fnames = joinpath.(data_dir,year_fnames[1:1])
ei,ej,ev,nodes = load_flowdata(fnames)
A = sparse(ei,ej,ev,length(nodes),length(nodes))

A = max.(A,A')


nnz(A)/size(A,1)
# fill!(A.nzval,1.0)
# d = vec(sum(A;dims=2))

using Plots
ENV["GKSwstype"] = "100"

# using StatsBase
# counts(A.nzval)

# m = minimum(A.nzval)
# cts = counts(A.nzval)
# plot(1:lastindex(cts),(1 .-cumsum(cts)./sum(cts))*(nnz(A)/size(A,1)),
#     xscale=:log10,leg=false)



#sample infrequent edges

#filter edges out 

#dont do anything 

B = deepcopy(A)
B = max.(B,B')
function downsample!(X::SparseMatrixCSC,val::Union{Int,Float64})
    X.nzval[X.nzval.<=val].=0
    dropzeros!(X)
end

B.nzval[B.nzval.<=20].=0
dropzeros!(B)

nnz(B)/size(B,1)

fill!(B.nzval,1.0)
B = max.(B,B')

# ncp,sets = make_ncp(B)

# ncp.size = map(c->min(c,size(A,1)-c),ncp.size)
# myncpplot1(ncp.size,ncp.cond)
# f = plot_hexbinncp1(ncp,A)

include("../code/graph-io.jl")
fnames

writeSMAT(B,joinpath(dst,"covidflows-2020_08_31-filtered-20.smat"))
fill!(A.nzval,1.0)

writeSMAT(A,joinpath(dst,"covidflows-2020_08_31.smat"))

