# figure for showing epidemic ncp sampling from a single diffusion

mainDir = "/p/mnt/scratch/network-epi/"
include(joinpath(mainDir,"code/fast-diffusion.jl"))
include(joinpath(mainDir,"code/graph-io.jl")) 
include(joinpath(mainDir,"code/data-io.jl")) 

using StatsPlots 
using Measures

#load data 
gname = "modmexico-city.smat"
A = loadGraph(gname,"input/graphs/")

#do a single epidemic 
E = EventEpidemic.SEIRData(A)

#single sample, many samples
Random.seed!(7)
seed_node = 2
zfigs = []
sfigs = []

z = zeros(Float64,size(A,1))
s = zeros(Float64,size(A,1))

l,E = EventEpidemic.epidemic(E,seed_node)

#single sample 
z.+=min.(E.itime,l+1)
s.+=E.snodes

pz = sweepcut(A,-z)
ps = sweepcut(A,-s)

f = StatsPlots.plot(pz.conductance,xscale=:log10,yscale=:log10,leg=false)
push!(zfigs,f)

f = StatsPlots.plot(ps.conductance,xscale=:log10,yscale=:log10,leg=false)
push!(sfigs,f)

for i = 1:20
    l,E = EventEpidemic.epidemic(E,seed_node)
    #stats
    z.+=min.(E.itime,l+1)
    s.+=E.snodes
end

pz = sweepcut(A,-z)
ps = sweepcut(A,-s)

f = StatsPlots.plot(pz.conductance,xscale=:log10,yscale=:log10,leg=false)
push!(zfigs,f)

f = StatsPlots.plot(ps.conductance,xscale=:log10,yscale=:log10,leg=false)
push!(sfigs,f)


@showprogress for i = 1:100
    l,E = EventEpidemic.epidemic(E,seed_node)
    z.+=min.(E.itime,l+1)
    s.+=E.snodes
end

pz = sweepcut(A,-z)
ps = sweepcut(A,-s)

f = StatsPlots.plot(pz.conductance,xscale=:log10,yscale=:log10,leg=false)
push!(zfigs,f)

f = StatsPlots.plot(ps.conductance,xscale=:log10,yscale=:log10,leg=false)
push!(sfigs,f)


#laying things out 
figs = vcat(sfigs,zfigs)
ybounds = vcat(ylims(figs[1])...)
for f in figs
    ybounds[1] = min(ylims(f)[1],ybounds[1])
    ybounds[2] = max(ylims(f)[2],ybounds[2])
end
for f in figs
    ylims!(f,(ybounds[1],ybounds[2]))
end

newf = Plots.plot(figs...,layout=(2,3),link=:both,size=(800,500))
Plots.plot!(newf[4],xlabel="Set Size",guidefontsize=15)
Plots.plot!(newf[4],ylabel="Conductance",guidefontsize=15,
        ylims=(ybounds[1],ybounds[2]),left_margin=4Measures.mm)
# Plots.plot!(newf,dpi=300)
# Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/example-figs/epidemic-ncp-sampling-comparison.png")



################################################################

#a better measure than conductance
# some ideas 
# conductance versus time (no degree normalization)
# (cut)
# (cut)/(inner set size)
# (cut)/(outer set size) - exapansion measure

function cutsize(A::SparseMatrixCSC,x::Vector)
    @assert(lastindex(A,1)==lastindex(A,2)==lastindex(x))
    #setting things up 
    cuts = zeros(Int,lastindex(A,1)-1)
    cut = 0

    colptr = A.colptr
    rowval = rowvals(A)
    nzval = A.nzval

    #sort nodes by ranking
    p = sortperm(x,rev=true) #p[1] node w/ highest ranking
    ranking = invperm(p)

    #use node ranking in place of a set S. 
    for (urank,u) in enumerate(p[1:end-1])
        for ind in colptr[u]:(colptr[u+1] - 1)
            v = rowval[ind]
            cut += ranking[v] <= urank ? -nzval[ind] : nzval[ind] 
        end
        cuts[urank] = cut
    end
    return cuts
end

function sweepexpansion(A::SparseMatrixCSC,x::Vector)
    @assert(lastindex(A,1)==lastindex(A,2)==lastindex(x))
    #setting things up #keep track of edges cut and node boundary
    setexpansions = zeros(Float64,lastindex(A,1)-1)
    cut = 0

    nodes = zeros(Bool,lastindex(A,1))
    nodecount = 0

    colptr = A.colptr
    rowval = rowvals(A)
    nzval = A.nzval

    #sort nodes by ranking
    p = sortperm(x,rev=true) #p[1] node w/ highest ranking
    ranking = invperm(p)

    #use node ranking in place of a set S. 
    for (urank,u) in enumerate(p[1:end-1])
        for ind in colptr[u]:(colptr[u+1] - 1)
            v = rowval[ind]

            if ranking[v] <= urank 
                cut -=  nzval[ind]
                nodecount -= nodes[v] 
                nodes[v] = 0
            else
                cut +=  nzval[ind]
                nodecount += 1-nodes[v]
                nodes[v] = 1
            end
            #remove node u from node boundary
            nodecount -= nodes[u] 
            nodes[u] = 0
            # @assert(sum(nodes)==nodecount)
        end
        setexpansions[urank] = cut/nodecount
    end
    return setexpansions
end

#########testing sweepexpansion
d = vec(sum(A;dims=2))
testoutput = sweepexpansion(A,d)

# #alternate computation 
ptest = sortperm(d,rev=true)
inds = ptest[1:227]
other = setdiff(1:lastindex(A,1),inds)

tmpcut = nnz(A[inds,other])
tmpnodeexpansion = count(x->x>0,vec(sum(A[inds,other];dims=1)))
count(x->x>0,sum(A[inds,other];dims=1))
testoutput[227] == tmpcut/tmpnodeexpansion


function sweepconductance(A::SparseMatrixCSC,x::Vector)
    @assert(lastindex(A,1)==lastindex(A,2)==lastindex(x))
    #setting things up 
    conds = zeros(Float64,lastindex(A,1)-1)
    cut = 0
    vol = 0 
    totalvol = nnz(A)

    colptr = A.colptr
    rowval = rowvals(A)
    nzval = A.nzval

    #sort nodes by ranking
    p = sortperm(x,rev=true) #p[1] node w/ highest ranking
    ranking = invperm(p)

    #use node ranking in place of a set S. 
    for (urank,u) in enumerate(p[1:end-1])
        for ind in colptr[u]:(colptr[u+1] - 1)
            v = rowval[ind]
            if ranking[v]<=urank
                cut -= nzval[ind]
            else
                cut+= nzval[ind]
            end
            vol += nzval[ind]
        end
        conds[urank] = cut/min(vol,totalvol-vol)
    end
    return conds
end

#########testing cutsize 
d = vec(sum(A;dims=2))
testoutput = cutsize(A,d)

# #computing another way 
ptest = sortperm(d,rev=true)
inds = ptest[1:227]
other = setdiff(1:lastindex(A,1),inds)
#final test conditon 
testoutput[227]==nnz(A[inds,other])

######testing sweepconductance
tvol = nnz(A)
testoutput = sweepconductance(A,d)
#alternate computation 
ptest = sortperm(d,rev=true)
inds = ptest[1:227]
other = setdiff(1:lastindex(A,1),inds)

tmpcut = nnz(A[inds,other])
tmpvol = nnz(A[inds,:])
tmpvol_other = nnz(A[other,:])

#final test condition
testoutput[227]==tmpcut/min(tmpvol,tvol-tmpvol_other)





## plot testing 
gname = "study-11-2023-0-noweak.smat"
gname = "study-11-2022-1.smat"
gname = "study-11-2022-50.smat"
gname = "modmexico-city.smat"
A = loadGraph(gname,"input/graphs/")


function single_plot(gname)
    @show gname
    A = loadGraph(gname,"pipeline/graphs/")
    E = EventEpidemic.SEIRData(A,beta=3e-1,gamma=1e-2)

    d = vec(sum(A;dims=2))
    seed_node = sortperm(d)[end-50]

    println("performing epidemics")
    z = zeros(Float64,size(A,1))
    l,E = EventEpidemic.epidemic(E,seed_node)
    #single sample 
    z.+=min.(E.itime,l+1)
    
    f = plot(sweepcut(A,-z).conductance,xscale=:log10,yscale=:log10,leg=false,ylims=(1e-4,1))
    f1 = plot(sweepcut(A,-z).cut,xscale=:log10,yscale=:log10,leg=false,ylims=(1,10^7))
    f2 = plot(sweepexpansion(A,-z),xscale=:log10,yscale=:log10,leg=false,ylims=(1e-2,1e2))

    @showprogress for i = 1:20
        l,E = EventEpidemic.epidemic(E,seed_node)
        #stats
        z.+=min.(E.itime,l+1)
    end

    plot!(f,sweepcut(A,-z).conductance,ylabel="conductance")
    plot!(f1,sweepcut(A,-z).cut,ylabel="cut")
    plot!(f2,sweepexpansion(A,-z),ylabel="expansion")

    return f,f1,f2
end


gnames = ["study-11-2023-0-noweak.smat",
    "study-11-2022-1.smat",
    # "study-11-2022-10.smat",
    # "study-11-2022-20.smat",
    # "study-11-2022-30.smat",
    # "study-11-2022-40.smat",
    "study-11-2022-50.smat"
]

figs = map(x->single_plot(x),gnames)
for i = 1:lastindex(gnames)
    Plots.show(figs[i][1])
end
plot!.(figs...,size=(300,200))


p1 = plot(first.(figs)...,layout=(1,7),size=(1500,200),ylabel="")
p2 = plot(map(x->x[2],figs)...,layout=(1,7),size=(1500,200),ylabel="")
p3 = plot(map(x->x[3],figs)...,layout=(1,7),size=(1500,200),ylabel="")

p1


p1 = plot(first.(figs)...,layout=(1,3),size=(800,200),ylabel="Conductance")
p2 = plot(map(x->x[2],figs)...,layout=(1,3),size=(800,200),ylabel="Cut")
p3 = plot(map(x->x[3],figs)...,layout=(1,3),size=(800,200),ylabel="Expansion")








res = single_plot(gname)
res[1]
E = EventEpidemic.SEIRData(A,beta=3e-1,gamma=5e-2)

z = zeros(Float64,size(A,1))
l,E = EventEpidemic.epidemic(E,seed_node)
#single sample 
z.+=min.(E.itime,l+1)
f = plot(sweepcut(A,-z).cut,xscale=:log10,yscale=:log10,leg=false,ylims=(1,nnz(A)))

@showprogress for i = 1:20
    l,E = EventEpidemic.epidemic(E,seed_node)
    #stats
    z.+=min.(E.itime,l+1)
end

plot!(f,sweepcut(A,-z).cut)






@time E = EventEpidemic.SEIRData(A,beta=8e-3,gamma=1e-5,tmax=200000);

l,E = EventEpidemic.epidemic(E,2)




#testing with study graphs 
gname = "study-11-2023-0-noweak.smat"
cap = 5000

A = loadGraph(gs[2],"input/graphs/")
E = EventEpidemic.SEIRData(A,beta=3e-2,qcap=cap);
res1 = []

@showprogress for i =1:50
    l,E = EventEpidemic.epidemic(E,rand(1:lastindex(A,1)),q=true);
    push!(res1,sum(EventEpidemic.get_infdata(l,E)[2]))
end
mean(res1)

A = loadGraph(gs[3],"input/graphs/")
E = EventEpidemic.SEIRData(A,beta=3e-2,qcap=cap);
res2 = []

@showprogress for i =1:50
    l,E = EventEpidemic.epidemic(E,rand(1:lastindex(A,1)),q=true);
    push!(res2,sum(EventEpidemic.get_infdata(l,E)[2]))
end
mean(res2)


A = loadGraph(gs[end],"input/graphs/")
E = EventEpidemic.SEIRData(A,beta=3e-2,qcap=cap);
res3 = []

@showprogress for i =1:50
    l,E = EventEpidemic.epidemic(E,rand(1:lastindex(A,1)),q=true);
    push!(res3,sum(EventEpidemic.get_infdata(l,E)[2]))
end
mean(res3)


scatter((vcat(mean(res1),mean(res2),mean(res3))./5e4).*100,leg=false)

# try qcap = (0,25,50,75,100,125,150,200,250,300,400,500,1000,2000,3000,4000,5000)





# x = MatrixNetworks.personalized_pagerank(A,0.9995,seed_node,1e-5)
# plot(MatrixNetworks.sweepcut(A,x).conductance,
#     yscale=:log10,xscale=:log10,leg=false)

# StatsPlots.plot(pz.conductance,xscale=:log10,yscale=:log10,leg=false)