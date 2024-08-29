#load in original graph and ncp rewired variant and stare at their ncps 


#script to perform ncp rewirings and save to scratch directory 
mainDir = joinpath("/",split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))]...)
include(joinpath(mainDir,"code/hierarchical-random-graphs/ncp-rewiring.jl"))
dstDir =  joinpath(mainDir,"scratch","ncp-rewired-graphs")
gpath = joinpath(mainDir,"pipeline/graphs/")

#gnames 
gs = [
        "study-20-draft-150",
        # "mexico-city",
        # "geometric",
        "study-25-1",
        "study-25-2",
        "study-25-150",
        "filtered",
        ]


gname = getgnames(gs[5],"input/graphs/")[1]
r_gname = "ncp-rewired-1.0-$gname"

g = canonical_graph_name(gname)[1:end-5]

A = loadGraph(gname,"pipeline/graphs/")
X = loadGraph(r_gname,"pipeline/graphs/")

ncp,_,sets = load_ncpdata(gname,"pipeline/data/$g/ncpdata/")
altncp,_,altsets = load_ncpdata(r_gname,"pipeline/data/$g/ncpdata/")

x = altncp.size
x = map(a->min(a,lastindex(A,1)-a),x)
y = altncp.cond

f = myhexbin_conditional(x,y,nbins=(80,60),
    ylims=(1e-3,1.2),
    xlims=(0.8,lastindex(A,1)))

myhexbin_conditional!(f,x,y,nbins=(80,60),
    ylims=(1e-3,1.2),
    xlims=(0.8,lastindex(A,1)))
