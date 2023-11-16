using ProgressMeter
using Distributed

#allows for execution from command line as well as an ide so files can be run modularlly
mainDir = joinpath("/",split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))]...)
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
