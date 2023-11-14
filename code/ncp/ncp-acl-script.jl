using ProgressMeter
using Distributed
# addprocs(25)

#allows for execution from command line as well as an ide so files can be run modularlly
mainDir = joinpath(split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))])
@everywhere mainDir = $mainDir

gpath = joinpath(mainDir,"input/graphs/")
@everywhere dstDir = joinpath(mainDir,"pipeline/data/")
#headless display mode for server when using GR backend for plotting. see the following
#https://discourse.julialang.org/t/generation-of-documentation-fails-qt-qpa-xcb-could-not-connect-to-display/60988
ENV["GKSwstype"] = "100"
@everywhere include(joinpath(mainDir,"code/graph-io.jl"))
@everywhere include(joinpath(mainDir,"code/ncp/ncp-acl.jl"))

gnames = readdir(gpath)
filter!(x->endswith(x,".smat"),gnames)

gnames = [getgnames(x,"pipeline/graphs/") for x in gnames]
for (i,x) in enumerate(gnames)
    if !startswith(x[1],"cn") && any(occursin.("cn-",x))
        filter!(x->!occursin("cn-",x),gnames[i])
    end
end

#sort gnames by number of edges 
p = sortperm(map(x->get_graph_stats(canonical_graph_name(x))[2],gnames))
gnames = gnames[p]

dsts = joinpath.(dstDir,map(x->first(x)[1:end-5],gnames),"ncpdata/")
dsts = repeat.([[x] for x in dsts],length.(gnames))

# testing that we got the right directories
# tests = first.(collect.(zip.(gnames,dsts)))
# x1 = map(x->x[1][1:end-5],tests)
# x2 = map(x->split(x[2],"/")[3],tests)
# all(x1.==x2)

ps = vcat(collect.(zip.(gnames,dsts))...) #[ (gname,gdst) ]

#filter out parameters where we've already done this computations
function check_graph(gname::String;dst::String="",all::Bool=true)
    if endswith(gname,".smat")
        g = gname[1:end-5]
    else
        g = gname[:]
    end

    fnames = readdir(dst)
    if all
        return ("ncpinfo-$g.txt" in fnames) && ("ncpmat-$g.smat" in fnames)
    else
        return ("ncpinfo-$g.txt" in fnames)    
    end
end
#keep parameters for graphs whose NCP we have not computed 
filter!(x->!check_graph(x[1],dst=x[2]),ps)
println("ngraphs for acl ncp: $(length(ps))")
#compute and save ncps
@showprogress pmap(x->make_ncp(x[1],gpath="pipeline/graphs/",dst=x[2],get_sets=true),ps)

#do flickr separately
# gnames = ["flickr"]
# gnames = [getgnames(x,"pipeline/graphs/") for x in gnames]
# for (i,x) in enumerate(gnames)
#     if !startswith(x[1],"cn") && any(occursin.("cn-",x))
#         filter!(x->!occursin("cn-",x),gnames[i])
#     end
# end
# dsts = joinpath.(dstDir,map(x->first(x)[1:end-5],gnames),"ncpdata/")
# dsts = repeat.([[x] for x in dsts],length.(gnames))
# ps = vcat(collect.(zip.(gnames,dsts))...) 


# @showprogress pmap(x->make_ncp(x[1],gpath="pipeline/graphs/",dst=x[2],get_sets=true),ps)


#custom calculation for study-20-draft-150 since the graph is small
# A = SparseMatrixCSC{Int,Int}(A)
# ncp_vec = []
# sets_vec = []
# @showprogress for i=1:1000  
#   ncp,sets = DiffusionAlgorithms.serial_weighted_ncp(A;
#           alpha=0.85, maxvisits=200, get_sets=true)
#   push!(ncp_vec,ncp)
#   push!(sets_vec,sets)
# end

# ncp = vcat(ncp_vec...)
# # sets = vcat(sets_vec...)
# using Measures
# f = myncpplot(map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1)),ncp.cond,nbins=60)
# Plots.plot!(f,xlabel="size",ylabel="conductance",xlims=(5,1.2*size(A,1)),
#     ylims=(1e-2,1),label="",margins=-5Measures.mm,
#     title="NCP - $(gname[1:end-5])\n ")
# title!(f,"")
# xlabel!(f,"")
# ylabel!(f,"")
# plot!(f,dpi=500)
# Plots.savefig(f,"ncp-$(gname[1:end-5]).png")


#generate all the plots 
# gnames = getgnames("cit-","input/graphs/")
# imgdsts = joinpath.(dstDir,map(x->first(x)[1:end-5],gnames),"imgs/")
# imgdsts = repeat.([[x] for x in imgdsts],length.(gnames))
# imgdsts = vcat(collect.(zip.(gnames,imgdsts))...)

# function load_ncp(gname,dataDir::String="data")
    
#     A = loadGraph(gname,"pipeline/graphs/")
#     n = nnz(A)

#     ncp,headerinfo = readdlm(joinpath(dataDir,"ncpinfo-$(gname[1:end-5]).txt"),',',header=true)
#     ncp = DataFrame(ncp,vec(headerinfo))

#     x = map(x->min(ncp.size[x],n-ncp.size[x]),1:size(ncp,1))
#     y = ncp.cond

#     return x,y
# end 

# include("ncpplots1.jl")


# #TODO test this out 
# @showprogress for (i,(gname,gdst)) in enumerate(ps)
#     #load data 
#     ncpsize,ncpcond = load_ncp(gname,gdst)
#     #make ncp
#     f = myncpplot1(ncpsize,ncpcond,nbins=100,xlims=(1,2e5),ylims=(1e-4,1.1))
#     plot!(f,title="ACL NCP\n$gname",xlabel="Size",ylabel="Conductance")
#     #save figure 
#     Plots.savefig(f,joinpath(imgdsts[i][2],"hexbinncp-$gname.png"))
# end


