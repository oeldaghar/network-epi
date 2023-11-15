#generating figure for appendix showing that removing subset of nodes is kosher

#basically, generate ncps for gname + compare w/ ncp for  modgname

mainDir = "/p/mnt/scratch/network-epi/code/"
include(joinpath(mainDir,"graph-io.jl"))
include(joinpath(mainDir,"data-io.jl"))

include(joinpath(mainDir,"ncp/ncp-acl.jl"))

using Plots
using Measures

gpath = "/homes/oeldagha/network-epidemics/all-graphs/graphs/"


function makegridlines!(f,xs,ys;color=:grey,linewidth=1,tickfontsize::Int=10)
  #jerry rig this to work for now. 
  #cant seem to find in the docs...

  #vertical lines
  for x in xs
      plot!(f,[10.0^x;10.0^x],[ylims(f)...];color,linewidth,leg=false,alpha=1)
      annotate!(f, (10.0^x-10.0^(x-0.85), ylims(f)[1], text("10^{$x}", tickfontsize, :right, RGB(0.0,0.0,0.0), :bottom)))
  end

  #horizontal lines
  for y in ys
      plot!(f,[xlims(f)...],[10.0^(y);10.0^(y)];color,linewidth,leg=false,alpha=1)
      annotate!(f, (1.1*xlims(f)[1], 10.0^y-10.0^(y-0.95),  text("10^{ $y}", tickfontsize, :left, RGB(0.0,0.0,0.0), :top)))
      # annotate!(f, (xlims(f)[1], 10.0^y,  text("10^{$y}", 9, :left, RGB(0.75,0.75,0.75), :top)))
  end

  #annotate axes
  # annotate!(f,(0.07*xlims(f)[2],1.8*ylims(f)[1],text(L"\textbf{Size}\rightarrow", 10, :left, RGB(0.25,0.25,0.25), :top)))
  # annotate!(f,(xlims(f)[1],20*ylims(f)[1],text(L"\leftarrow\textbf{Conductance}", 20, :left, RGB(0.25,0.25,0.25), :top)))

  return f
end


function make_acl_plot(A,ncp) 
  x,y = ncp.size,ncp.cond
  nnodes = size(A,1)

  x = map(a->min(a,nnodes-a),x)
  xbounds,ybounds = (1.0,2e5),(2e-6,1.15)

  f = plot(margin=-20Measures.mm,  topmargin=isempty("") ? 0Measures.mm : -nlines*0.5Measures.mm , size=(600,400), 
                  ylims=ybounds, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
  makegridlines!(f, [1,2,3,4,5], -(1:5),tickfontsize=15)

  c = cgrad(cgrad(:viridis)[1:250])

  myhexbin2!(f,x,y,xlims=xbounds,ylims=ybounds,color=c,nbins=80,normalize_color=true)
  myhexbin2!(f,x,y,xlims=xbounds,ylims=ybounds,color=c,nbins=80,normalize_color=true)
  return f 
end

#mexico city
# gname = "mexico-city.smat"
# A = loadGraph(gname,gpath)
# A = SparseMatrixCSC{Int,Int}(A)

# println("working on ncp for $gname")
# dt = @elapsed ncp,sets = DiffusionAlgorithms.serial_weighted_ncp(A;
#         alpha=0.995, maxvisits=30, get_sets=false, maxlarge = 15,
#         epsbig=1e-5);

# println("ncp completed for $gname\nelapsed time: $dt seconds")

# #save data
# dst = "/p/mnt/scratch/network-epi/pipeline/deep-ncp-data/"
# # CSV.write(joinpath(dst,"ncpinfo-$(gname[1:end-5]).txt"),ncp)

# #load from file 
# ncp,headerinfo = readdlm(dst*"ncpinfo-$(gname[1:end-5]).txt",',',header=true)
# ncp = DataFrame(ncp,vec(headerinfo))

# x,y = ncp.size,ncp.cond
# nnodes = size(A,1)

# x = map(a->min(a,nnodes-a),x)
# xbounds,ybounds = (1.0,2e5),(2e-6,1.15)

# f = plot(margin=-20Measures.mm,  topmargin=isempty("") ? 0Measures.mm : -nlines*0.5Measures.mm , size=(600,400), 
#                 ylims=ybounds, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
# makegridlines!(f, [1,2,3,4,5], -(1:5),tickfontsize=15)

# c = cgrad(cgrad(:viridis)[1:250])

# myhexbin2!(f,x,y,xlims=xbounds,ylims=ybounds,color=c,nbins=80,normalize_color=true)
# myhexbin2!(f,x,y,xlims=xbounds,ylims=ybounds,color=c,nbins=80,normalize_color=true)

# plot!(f,colorbar=true,
#         bottom_margin=-5Measures.mm,
#         right_margin=-1Measures.mm,
#         dpi=1000)
# savefig(f,"/p/mnt/scratch/network-epi/code/paper-figs/large-set-removal/acl-ncp-$(gname[1:end-5]).png")


# #mod mexico city on same axis for reference
# gname = "modmexico-city.smat"
# A = loadGraph("$gname","input/graphs/")
# dloc = "/p/mnt/scratch/network-epi/pipeline/data/$(gname[1:end-5])/ncpdata/"
# ncp,headerinfo = readdlm(dloc*"ncpinfo-$(gname[1:end-5]).txt",',',header=true)
# ncp = DataFrame(ncp,vec(headerinfo))

# x,y = ncp.size,ncp.cond
# nnodes = size(A,1)

# x = map(a->min(a,nnodes-a),x)
# xbounds,ybounds = (1.0,2e5),(2e-6,1.15)

# f = plot(margin=-20Measures.mm,  topmargin=isempty("") ? 0Measures.mm : -nlines*0.5Measures.mm , size=(600,400), 
#                 ylims=ybounds, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
# makegridlines!(f, [1,2,3,4,5], -(1:5),tickfontsize=15)

# c = cgrad(cgrad(:viridis)[1:250])

# myhexbin2!(f,x,y,xlims=xbounds,ylims=ybounds,color=c,nbins=80,normalize_color=true)
# myhexbin2!(f,x,y,xlims=xbounds,ylims=ybounds,color=c,nbins=80,normalize_color=true)

# plot!(f,colorbar=true,
#         bottom_margin=-5Measures.mm,
#         right_margin=-1Measures.mm,
#         dpi=1000)
# savefig(f,"/p/mnt/scratch/network-epi/code/paper-figs/large-set-removal/acl-ncp-$(gname[1:end-5]).png")



#=

#fb networks 
gpath = "/p/mnt/data/raw-data/porter-facebook/facebook100/"
# modUIllinois20
gname = "UIllinois20.smat"
A = loadGraph(gname,gpath)
A = SparseMatrixCSC{Int,Int}(A)

println("working on ncp for $gname")
dt = @elapsed ncp,sets = DiffusionAlgorithms.serial_weighted_ncp(A;
        alpha=0.995, maxvisits=30, get_sets=false, maxlarge = 15,
        epsbig=1e-5);

println("ncp completed for $gname\nelapsed time: $dt seconds")


#save data
dst = "/p/mnt/scratch/network-epi/pipeline/deep-ncp-data/"
CSV.write(joinpath(dst,"ncpinfo-$(gname[1:end-5]).txt"),ncp)


# modWisconsin87
gname = "Wisconsin87.smat"
A = loadGraph(gname,gpath)
A = SparseMatrixCSC{Int,Int}(A)

println("working on ncp for $gname")
dt = @elapsed ncp,sets = DiffusionAlgorithms.serial_weighted_ncp(A;
        alpha=0.995, maxvisits=30, get_sets=false, maxlarge = 15,
        epsbig=1e-5);
println("ncp completed for $gname\nelapsed time: $dt seconds")

#save data
dst = "/p/mnt/scratch/network-epi/pipeline/deep-ncp-data/"
CSV.write(joinpath(dst,"ncpinfo-$(gname[1:end-5]).txt"),ncp)
=#


### deep dive in study graphs ###
#= 
gpath = "input/graphs/"
gs = [ "study-11-2023-0-noweak.smat",
    "study-11-2022-1.smat",
    "study-11-2022-10.smat",
    "study-11-2022-20.smat",
    "study-11-2022-30.smat",
    "study-11-2022-40.smat",
    "study-11-2022-45.smat",
    "study-11-2022-50.smat",
]

for gname in gs 
        A = loadGraph(gname,gpath)
        A = SparseMatrixCSC{Int,Int}(A)

        println("working on ncp for $gname")
        dt = @elapsed ncp,sets = DiffusionAlgorithms.serial_weighted_ncp(A;
                alpha=0.95, maxvisits=50, get_sets=false, maxlarge = 50,
                epsbig=1e-10);

        println("ncp completed for $gname\nelapsed time: $dt seconds")

        dst = "/p/mnt/scratch/network-epi/pipeline/deep-ncp-data/"
        CSV.write(joinpath(dst,"ncpinfo-$(gname[1:end-5]).txt"),ncp)
end

dst = "/p/mnt/scratch/network-epi/pipeline/deep-ncp-data/"

#looking at results 
gname = "mexico-city.smat"
ncp,headerinfo = readdlm(dst*"ncpinfo-$(gname[1:end-5]).txt",',',header=true)
ncp = DataFrame(ncp,vec(headerinfo))

x,y = ncp.size,ncp.cut
nnodes = 5e4

x = map(a->min(a,nnodes-a),x)
xbounds,ybounds = (1.0,2e5),(2e-6,1.15)

f = plot(margin=-20Measures.mm,  topmargin=isempty("") ? 0Measures.mm : -nlines*0.5Measures.mm , size=(600,400), 
                ylims=ybounds, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
makegridlines!(f, [1,2,3,4,5], -(1:5),tickfontsize=15)

c = cgrad(cgrad(:viridis)[1:250])

myhexbin2!(f,x,y,xlims=xbounds,ylims=ybounds,color=c,nbins=80,normalize_color=true)
myhexbin2!(f,x,y,xlims=xbounds,ylims=ybounds,color=c,nbins=80,normalize_color=true)



gnames = getgnames("study-weak","input/graphs/")
gname = "study-noweak.smat"

gnames[5]

A = loadGraph(gnames[5],"input/graphs/")
A = SparseMatrixCSC{Int,Int}(A)

println("working on ncp for $gname")
dt = @elapsed ncp,sets = DiffusionAlgorithms.serial_weighted_ncp(A;
        alpha=0.95, maxvisits=20, get_sets=false, maxlarge = 10,
        epsbig=1e-4);


x,y = ncp.size,ncp.cond
nnodes = lastindex(A,1)

x = map(a->min(a,nnodes-a),x)
xbounds,ybounds = (1.0,2e5),(2e-6,1.15)

f = plot(margin=-20Measures.mm,  topmargin=isempty("") ? 0Measures.mm : -nlines*0.5Measures.mm , size=(600,400), 
                ylims=ybounds, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
makegridlines!(f, [1,2,3,4,5], -(1:5),tickfontsize=15)

c = cgrad(cgrad(:viridis)[1:250])

myhexbin2!(f,x,y,xlims=xbounds,ylims=ybounds,color=c,nbins=80,normalize_color=true)
myhexbin2!(f,x,y,xlims=xbounds,ylims=ybounds,color=c,nbins=80,normalize_color=true)

=#
xbounds,ybounds = (1.0,2e5),(2e-6,1.15)
c = cgrad(cgrad(:viridis)[1:250])

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

#do 100 ppr and plot the dist 

# gnames = getgnames("long","input/graphs/")

gnames = ["study-11-2023-0-noweak.smat",
"study-11-2022-1.smat",
"study-11-2023-1-longrange-1.smat",
"study-11-2023-1-longrange-2.smat",
"study-11-2023-1-longrange-3.smat",
"study-11-2023-1-longrange-5.smat",
"study-11-2023-1-longrange-8.smat"]


gnames = filter!(x->endswith(x,".smat"),readdir("input/tmp/"))[1:end-1]
sort!(gnames,by=x->parse(Int,split(x[1:end-5],"-")[end]))

gnames = gnames[[1;2;11;21]]
gnames = getgnames("test","input/tmp/")
figs = []
@showprogress for gname in gnames
        A = loadGraph(gname,"input/tmp/")
        A = SparseMatrixCSC{Int,Int}(A)

        xedges = vcat(0:1:99,100:10:1000-1,1000:50:1e4-1,1e4:5e2:1e5)
        yedges = 10.0 .^(range(-6,0,120))
        dvec = vec(sum(A;dims=2))
        using StatsBase
        #initialize histogram 
        h = fit(Histogram,([0],[0]),(xedges,yedges))

        @showprogress for i=1:100
                ppr_x = MatrixNetworks.personalized_pagerank(A,0.995,rand(1:lastindex(A,1)),1e-2);
                p = MatrixNetworks.sweepcut(A,ppr_x./dvec);
                StatsBase.merge!(h,fit(Histogram,(1:lastindex(p.conductance),p.conductance),(xedges,yedges)))
        end

        xres,yres = histogram_to_data(h)

        f = plot()
        # myhexbin_conditional!(f,xres,yres,xlims=xbounds,ylims=(1e-2,2),color=c,nbins=80)
        myhexbin2!(f,xres,yres,xlims=xbounds,ylims=(1e-2,2),color=c,nbins=80,normalize_color=true)
        myhexbin2!(f,xres,yres,xlims=xbounds,ylims=(1e-2,2),color=c,nbins=80,normalize_color=true)
        plot!(f,xticks=[1,10,100,1000,1e4,1e5],title=gname[1:end-5])
        push!(figs,f)
end

B = loadGraph(gnames[1],"input/tmp/")
nnz(B)/lastindex(B,1)
figs[1]



figs[2]
figs[3]
figs[4]

figs[5]
figs[6]
figs[7]
figs[8]
figs[9]
figs[10]
figs[11]
figs[21]

A = loadGraph(gnames[1],"input/tmp/")
B = loadGraph(gnames[2],"input/tmp/")
norm(A.-B)

##############################
include("../fast-diffusion.jl")
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

gnames = ["study-11-2023-0-noweak.smat",
"study-11-2022-1.smat",
"study-11-2023-1-longrange-1.smat",
"study-11-2023-1-longrange-2.smat",
"study-11-2023-1-longrange-3.smat",
"study-11-2023-1-longrange-5.smat",
"study-11-2023-1-longrange-8.smat",
"study-11-2023-1-longrange-10.smat",
"study-11-2023-1-longrange-12.smat",
"study-11-2023-1-longrange-15.smat",]


figs = []
@showprogress for gname in gnames
        A = loadGraph(gname,"input/tmp/")
        A = SparseMatrixCSC{Int,Int}(A)

        E = EventEpidemic.SEIRData(A,beta=1e-1)

        xedges = vcat(0:1:99,100:10:1000-1,1000:50:1e4-1,1e4:5e2:1e5)
        yedges = 10.0 .^(range(-6,0,120))
        dvec = vec(sum(A;dims=2))
        using StatsBase
        #initialize histogram 
        h = fit(Histogram,([0],[0]),(xedges,yedges))

        @showprogress for i=1:100
                xrank = epi_sample(E,3)
                p = MatrixNetworks.sweepcut(A,xrank);
                StatsBase.merge!(h,fit(Histogram,(1:lastindex(p.conductance),p.conductance),(xedges,yedges)))
        end

        #get data for hexbin plotting 
        xres,yres = histogram_to_data(h)

        f = plot()
        # myhexbin_conditional!(f,xres,yres,xlims=xbounds,ylims=(1e-2,2),color=c,nbins=80)
        myhexbin2!(f,xres,yres,xlims=xbounds,ylims=(1e-2,2),color=c,nbins=80,normalize_color=true)
        myhexbin2!(f,xres,yres,xlims=xbounds,ylims=(1e-2,2),color=c,nbins=80,normalize_color=true)
        plot!(f,xticks=[1,10,100,1000,1e4,1e5],title=gname[1:end-5])
        push!(figs,f)
end

figs[1]
figs[2]
figs[3]
figs[4]
figs[5]
figs[6]
figs[7]
figs[8]
figs[9]
figs[10]



##testing w/ sbm 
import Graphs: stochastic_block_model
using Random

# p1 = 20
# p2 = 5e-4
# n = Int(1e4)
# Random.seed!(27)
# A = sparse(stochastic_block_model([p1 p2 0; p2 p1 p2; 0 p2 p1],[n, n, n]))

# A = sparse(erdos_renyi_undirected(3*n,20))


# A = loadGraph("modmexico-city.smat","input/graphs/")
A = loadGraph("study-11-2023-1-longrange-15.smat","input/graphs/")

include("../fast-diffusion.jl")

E = EventEpidemic.SEIRData(A,beta=9e-3)

xedges = vcat(0:1:99,100:10:1000-1,1000:50:1e4-1,1e4:5e2:1e5)
yedges = 10.0 .^(range(-6,0,120))
dvec = vec(sum(A;dims=2))
using StatsBase
#initialize histogram 
h = fit(Histogram,([0],[0]),(xedges,yedges))

@showprogress for i=1:100
        xrank = epi_sample(E,3)
        p = MatrixNetworks.sweepcut(A,xrank);
        StatsBase.merge!(h,fit(Histogram,(1:lastindex(p.conductance),p.conductance),(xedges,yedges)))
end

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

f = plot()
# myhexbin_conditional!(f,xres,yres,xlims=xbounds,ylims=(1e-2,2),color=c,nbins=80)
myhexbin2!(f,xres,yres,xlims=xbounds,ylims=(1e-2,2),color=c,nbins=80,normalize_color=true)
myhexbin2!(f,xres,yres,xlims=xbounds,ylims=(1e-2,2),color=c,nbins=80,normalize_color=true)
plot!(f,xticks=[1,10,100,1000,1e4,1e5],title="SBM Model")




minx,miny = get_approximate_mins(xres,yres,50000)
nnodes = lastindex(A,1)
extrapolate!(minx,miny,2*nnodes)
#compute auc from y=1 to conductance minimum
# auc = approximate_auc(log10.(minx), log10.(1 .-miny)) #auc from y=1 to min conductance
auc = approximate_auc(log10.(minx), -1 .*log10.(1 .-miny))
auc/log10(nnodes) #normalizing 


1-approximate_auc(minx./nnodes,miny)



seed_node = rand(1:lastindex(A,1))
l,E = EventEpidemic.epidemic(E,seed_node)
using DelimitedFiles
Axy = readdlm("input/graphs/study-11-2023-1-longrange-10.xy",'\t')
Axy


scatter(Axy[:,1],Axy[:,2],markersize=0.5,markerstrokewidth=0,leg=false,size=(900,600))