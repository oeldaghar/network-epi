#allows for execution from command line as well as an ide so files can be run modularlly
mainDir = joinpath("/",split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))]...)

include(joinpath(mainDir,"code/graph-io.jl")) 
include(joinpath(mainDir,"code/ncp/ncpplots1.jl")) 
include(joinpath(mainDir,"code/ncp/ncp-acl.jl"))

using DelimitedFiles
using DataFrames
using Measures
using Pkg
using PerceptualColourMaps
using LaTeXStrings
using ProgressMeter
##
gr()
# ENV["GRDIR"] = ""
# Pkg.build("GR")
#functions for loading in data 

function graphSize(gname)
    A = loadGraph(joinpath("pipeline/graphs/",gname))
    return size(A,1)
end

##
function load_ncp(gname,datadir="pipeline/data/",normalizex::Bool=true;ncptype="acl")
    
    base_graph = canonical_graph_name(gname)
    ncptype = lowercase(ncptype)

    #load info 
    if ncptype == "acl"
        ncp,headerinfo = readdlm(joinpath(datadir,base_graph[1:end-5],"ncpdata","ncpinfo-$(gname[1:end-5]).txt"),',',header=true)
    elseif ncptype=="epi" #epidemic ncp
        ncp,headerinfo = readdlm(joinpath(datadir,base_graph[1:end-5],"ncpdata","ncpinfo-epidemic-subsampling-4-seir-50000-$(gname[1:end-5]).txt"),',',header=true)
    else
        AssertionError("ncptype must be either: `acl' or `epi'")
    end
    ncp = DataFrame(ncp,vec(headerinfo))

    #use smaller sized set.
    if normalizex 
        n = graphSize(gname)
        x = map(x->min(ncp.size[x],n-ncp.size[x]),1:size(ncp,1))
    else 
        x = ncp.size
    end

    y = ncp.cond

    return x,y
end 


function makegridlines!(f,xs,ys;color=:grey,linewidth=1,tickfontsize::Int=10)
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
    return f
end

function make_base_plot(nlines,xbounds=(0.8,2e5),tickfontsize=15)
    ydim = (nlines+1)*35 + (nlines == 2 ? 1 : 0)
    ylims = (10.0^(-nlines-1),2.0)
    @show ydim, ylims
    #perform plotting 
    f = plot(margin=-20Measures.mm,  top_margin= 0*Measures.mm, size=(250,ydim), 
              ylims=ylims, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)

    makegridlines!(f, [1,2,3,4,5], -(1:nlines),tickfontsize=tickfontsize)

    return f 
end

#epidemic ncp plot
function make_epi_fig(gname; nlines = 3, xbounds = (0.8,2e5), nbins=(80,60),
                  color=:inferno,clims=(0,1),tickfontsize=15)
    
    println(gname)
    f = make_base_plot(nlines,xbounds,tickfontsize)

    x,y = load_ncp(gname,ncptype="epi")

    ybounds = (10.0^(-nlines-1),2.0)
    #place ncp data on plot 
    myhexbin_conditional!(f,x,y;
        xlims=xbounds,ylims=ybounds,nbins=nbins,
        color=color)
    #plotting a second time 
    myhexbin_conditional!(f,x,y;
        xlims=xbounds,ylims=ybounds,nbins=nbins,
        color=color)
    
    plot!(f;ylims=ybounds,clims=clims)

    return f
end


#stand alone acl ncp plots
function make_acl_fig(gname;
                    nlines = 3, 
                    # ybounds = (1e-3,1),
                    xbounds = (0.8,2e5),
                    tickfontsize=15,
                    nbins=60,
                    color=cgrad(cgrad(:viridis)[1:250]))
    
    
    f = make_base_plot(nlines,xbounds,tickfontsize)
    ybounds = (10.0^(-nlines-1),2.0)

    x,y = load_ncp(gname,ncptype="acl")
    myhexbin2!(f,x,y,
            # xlims=xbounds,ylims=ybounds,
            color=color,nbins=nbins,normalize_color=true)
    myhexbin2!(f,x,y,
            # xlims=xbounds,ylims=ybounds,
            color=color,nbins=nbins,normalize_color=true)
    
    plot!(f;ylims=ybounds,xlims=xbounds,clims=(0,1))
    return f
end


# missed sets functions 
function _plot_missed_sets1!(f,ncp::DataFrame,ncpmat::SparseMatrixCSC,snodes,sets;
            nbins::Int=100,ntrials::Int=1,color=:inferno,xbounds::Tuple=(1,2e5),ybounds::Tuple=(3e-5,1))
   
    @assert(maximum(snodes)<=ntrials,"ntrials < max # of times a node is susceptible")
    #larger values means we missed more nodes in that set
    missedsets = ncpmat*snodes./sum(ncpmat;dims=2)
    #for aggregating over different diffusions
    missedsets./=ntrials

    ncp.size = map(x->min(x,lastindex(snodes)-x),ncp.size)
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
        # nodes = collect(union(sets[s]...))
        # bin_val = (sum(snodes[nodes])/length(nodes))/ntrials 
        #writing this as a generator instead 
        
        node_vec = zeros(Bool,lastindex(snodes))
        for sind in s
            tmpset = sets[sind]
            for node in tmpset
                node_vec[node] = 1
            end
        end

        bin_val = (sum(snodes[node_vec])/sum(node_vec))/ntrials 

        for k=1:8
            push!(wh,bin_val)
        end
    end

    Plots.plot!(f,10.0.^xsh,10.0.^ysh,fill_z=wh,linecolor=nothing,
        seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=false,color=color)
    
    Plots.plot!(f,10.0.^xsh,10.0.^ysh,fill_z=wh,linecolor=nothing,
        seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=false,color=color)
    
    return f
end

function _plot_missed_sets1!(f,gname::String;ncploc::String="",sloc::String="",
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
        f_loc = _plot_missed_sets1!(f_loc,ncp,ncpmat,snodes_col,sets;nbins=nbins,ntrials=ntrials,color=color)
        Plots.plot!(f_loc,title="$gname - $model($beta,0.05,$(qpercent) Qpercent)\nAvg. Infected Fraction: $avg_inf_frac%",colorbar_title="Fraction Nodes Susceptible")
        figs[(beta,qpercent)] = deepcopy(f_loc) #(beta,qpercent): missed sets fig
    end
    return figs 
end


function make_missed_sets_figs(gname; nlines = 3, 
                    # ybounds = (1e-3,1),
                    xbounds = (1.0,2e5),
                    tickfontsize=15,
                    nbins=60,
                    color=:inferno,
                    betas::Vector{Float64} = collect(0.02:0.02:0.1))
    
    println(gname)

    #path information and naming
    g = canonical_graph_name(gname)
    ncploc = "pipeline/data/$(g[1:end-5])/ncpdata/"
    sloc = "pipeline/data/$(g[1:end-5])/diffusions/uniform/"
    
    #set down gridlines 
    f = make_base_plot(nlines,xbounds,tickfontsize)
    
    #missed sets plots    
    tmpfs = _plot_missed_sets1!(deepcopy(f),gname,ncploc=ncploc,sloc=sloc,ntrials=50,model="seir",
                betas=betas,nbins=nbins,color=color)

    return tmpfs
end


function makefigs(gname,betas,nlines::Int=3,nbins=60)
    xbounds = (0.8,2e5)
    tickfontsize = 10
    #specify colors 
    missed_sets_c = cgrad(cgrad(:magma)[1:240])
    acl_c = cgrad(cgrad(:viridis)[1:250])
    epi_c = cgrad(cmap("CBL2")[75:230])


    #missed sets 
    tmpfigs = make_missed_sets_figs(gname,nlines=nlines,betas=betas,
                        color=missed_sets_c,
                        nbins=nbins,
                        xbounds=xbounds,
                        tickfontsize=tickfontsize)
    fs = [tmpfigs[(beta,qpercent)] for (beta,qpercent) in sort(collect(keys(tmpfigs)))]
    plot!.(fs,size=(600,400),clims=(0,1),title="")

    #acl 
    acl_f = make_acl_fig(gname,nlines=nlines,nbins=nbins,
                        xbounds=xbounds,
                        color=acl_c,
                        tickfontsize=tickfontsize)
    plot!(acl_f,size=(600,400),clims=(0,1))

    #epidemic
    epi_f = make_epi_fig(gname,nlines=nlines,xbounds=xbounds,color=epi_c,
            tickfontsize=tickfontsize)
    plot!(epi_f,size=(600,400),clims=(0,1))


    return vcat([epi_f],[acl_f],fs)
end



### specifying parameters 
#gname,betas,nlines
hs = ["commutes-all", "mexico","filtered", #row1
    "cn-moduillinois", "cn-Penn", "cn-modWiscon", #row2...
    "dblp", "enron", "anon", #row 3
    "cit-HepPh", "slashdot", #"flickr", #skip flickr for now 
    "geometric","study-25-150","cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.9-0.0-100-5.smat"]
    #  "rewired-10000.0-modmexico","rewired-10000.0-cn-moduillinois","er-10000.0-dblp-cc"]
gnames = map(h->h=="" ? "" : getgnames(h,"pipeline/graphs/")[1], hs)

beta_vec = [
    vcat(1e-3:1e-3:5e-3), #commutes 
    vcat(2e-2:2e-2:1e-1), #mexico 
    [0.004, 0.006, 0.008, 0.01, 0.03], #covid flows 
    [0.04, 0.07, 0.1, 0.13, 0.16], #uill 
    [0.04, 0.07, 0.1, 0.13, 0.16], #penn
    [0.04, 0.07, 0.1, 0.13, 0.16], #wisc
    vcat(2e-2:2e-2:1e-1), #dblp
    [0.006, 0.008, 0.01, 0.03, 0.05], #enron
    vcat(2e-2:2e-2:1e-1), #anon 
    [0.004, 0.006, 0.008, 0.01, 0.03], #cit 
    [0.003, 0.005, 0.007, 0.009, 0.02], #slashdot 
    vcat(7e-3:1e-3:1e-2,2e-2), #geometric 
    [2e-3, 4e-3, 5e-3, 6e-3, 8e-3], #geo-comm 
    [0.009, 0.02, 0.04, 0.06, 0.08], #rwc 
]



#triming the middle beta
beta_vec = [vcat(x[1:2],x[4:5]) for x in beta_vec]


nlines_vec = [
    2, #commutes 
    4, #mex 
    3, #flows 
    2, 2, 2, #sparse fb 
    3, #dblp 
    3, #email
    3, #anon
    2, #cit
    2, #slashdot 
    3, #geometric
    2, #geo-comm  
    4, #rwc 
]


##testing 
# ind = 2
# gname = gnames[ind]
# # test_fs
# test_fs = makefigs(gnames[ind],beta_vec[ind],nlines_vec[ind],40)

# plot!.(test_fs,size=(300,200),ylims=(ylims(test_fs[1])[1],1.5))
# newf = plot(test_fs...,layout=(1,lastindex(test_fs)),size=(300*lastindex(test_fs),200),
#     margins=0Measures.mm,
#     left_margin=-12Measures.mm,
#     right_margin=-4Measures.mm,
#     top_margin=-2Measures.mm,
#     bottom_margin=-3Measures.mm)
# plot!(newf[1],right_margin=0.5Measures.mm)
# plot!(newf[2],right_margin=0.5Measures.mm)


for ind = 1:lastindex(gnames)
    gname = gnames[ind]
    test_fs = makefigs(gnames[ind],beta_vec[ind],nlines_vec[ind],40)

    plot!.(test_fs,size=(300,200),ylims=(ylims(test_fs[1])[1],1.5))
    newf = plot(test_fs...,layout=(1,lastindex(test_fs)),size=(300*lastindex(test_fs),200),
        margins=0Measures.mm,
        left_margin=-12Measures.mm,
        right_margin=-4Measures.mm,
        top_margin=-2Measures.mm,
        bottom_margin=-3Measures.mm)
    plot!(newf[1],right_margin=0.5Measures.mm)
    plot!(newf[2],right_margin=0.5Measures.mm)
    
    plot!(newf,dpi=500)
    Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets-v2-$(gname[1:end-5]).png")
end





gname = getgnames("study-25-150","input/graphs/")[1]
test_fs = makefigs(gname,[0.02],4,80)
plot!.(test_fs,size=(300,200),ylims=(ylims(test_fs[1])[1],1.5))

newf = plot(test_fs...,layout=(1,lastindex(test_fs)),size=(300*lastindex(test_fs),200),
    margins=0Measures.mm,
    left_margin=0Measures.mm,
    right_margin=0Measures.mm,
    top_margin=0Measures.mm,
    bottom_margin=0Measures.mm)


gnames = getgnames("study-25-150","input/graphs/")
# gnames = getgnames("11-2023-0","input/graphs/")
acl_c = cgrad(cgrad(:viridis)[1:250])
acl_f = make_acl_fig(gnames[1],nlines=4,nbins=50,
    xbounds=xbounds,
    color=acl_c,
    tickfontsize=15)
plot!(acl_f,size=(600,400),clims=(0,1))
title!("study-150")