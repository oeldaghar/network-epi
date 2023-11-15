using Distributed
gpath = "input/graphs/"
dstDir = "pipeline/data/"
parentDir = "/p/mnt/scratch/network-epi/"


using PerceptualColourMaps #for colorblind colormaps 
using Plots

using Measures
#headless display mode for server when using GR backend for plotting. see the following
#https://discourse.julialang.org/t/generation-of-documentation-fails-qt-qpa-xcb-could-not-connect-to-display/60988
ENV["GKSwstype"] = "100"

@everywhere include(joinpath(parentDir,"code/graph-io.jl"))
@everywhere include(joinpath(parentDir,"code/data-io.jl"))
@everywhere include(joinpath(parentDir,"code/ncp/ncp-acl.jl"))
# gnames = readdir(gpath)
# filter!(x->endswith(x,".smat"),gnames)

# gnames = [getgnames(x,"pipeline/graphs/") for x in gnames]
# for (i,x) in enumerate(gnames)
#     if !startswith(x[1],"cn") && any(occursin.("cn-",x))
#         filter!(x->!occursin("cn-",x),gnames[i])
#     end
#     filter!(x->!occursin("triangle-rewired-",x),gnames[i])
# end



# dataDir = joinpath.(dstDir,map(x->first(x)[1:end-5],gnames))

# dsts = joinpath.(dstDir,map(x->first(x)[1:end-5],gnames),"ncpdata/")
# imgdsts = joinpath.(dstDir,map(x->first(x)[1:end-5],gnames),"imgs/missed-sets/")
# slocs = joinpath.(dstDir,map(x->first(x)[1:end-5],gnames),"diffusions/uniform/")

# dsts = repeat.([[x] for x in dsts],length.(gnames))
# imgdsts = repeat.([[x] for x in imgdsts],length.(gnames))
# slocs = repeat.([[x] for x in slocs],length.(gnames))

# pathinfo = vcat(collect.(zip.(gnames,dsts,slocs,imgdsts))...) 



# #filter out parameters where we've already done this computations
# function check_graph(gname::String;dst::String="",all::Bool=true)
#     if endswith(gname,".smat")
#         g = gname[1:end-5]
#     else
#         g = gname[:]
#     end

#     fnames = readdir(dst)
#     if all
#         return ("ncpinfo-$g.txt" in fnames) && ("ncpmat-$g.smat" in fnames)
#     else
#         return ("ncpinfo-$g.txt" in fnames)    
#     end
# end
# #keep parameters for graphs whose NCP we have not computed 
# filter!(x->check_graph(x[1],dst=x[2]),pathinfo)
# println("ngraphs for acl ncp: $(length(pathinfo))")

# @showprogress for (gname,ncploc,sloc,imgdst) in pathinfo
#     for model in ["seir"]
#         println("working on $model for $gname")
#         figs = plot_missed_sets1(gname,ncploc=ncploc,sloc=sloc,ntrials=50,model=model)
        
#         for (beta,qpercent) in keys(figs)
#             Plots.savefig(figs[(beta,qpercent)],joinpath(imgdst,"$(gname[1:end-5])-$model-$beta-0.05-qpercent-$qpercent-missed_sets1.png"))
#         end
#     end
# end


####################################################################################################
############################################################################################
#####################################################################

function graphSize(gname)
    A = loadGraph(joinpath("pipeline/graphs/",gname))
    return size(A,1)
end

# function makegridlines!(f,xs,ys;color=:grey,linewidth=1,tickfontsize::Int=10)
#     #jerry rig this to work for now. 
#     #cant seem to find in the docs...

#     #vertical lines
#     for x in xs
#         plot!(f,[10.0^x;10.0^x],[ylims(f)...];color,linewidth,leg=false,alpha=1)
#         plot!(f,[10.0^x;10.0^x],[ylims(f)...];color,linewidth,leg=false,alpha=1)
#         annotate!(f, (10.0^x-10.0^(x-0.85), ylims(f)[1], text("10^{$x}", tickfontsize, :right, RGB(0.0,0.0,0.0), :bottom)))
#     end

#     #horizontal lines
#     for y in ys
#         plot!(f,[xlims(f)...],[10.0^(y);10.0^(y)];color,linewidth,leg=false,alpha=1)
#         plot!(f,[xlims(f)...],[10.0^(y);10.0^(y)];color,linewidth,leg=false,alpha=1)
#         annotate!(f, (1.1*xlims(f)[1], 10.0^y+10^(y+0.2),  text("10^{ $y}", tickfontsize, :left, RGB(0.0,0.0,0.0), :top)))
#         # annotate!(f, (xlims(f)[1], 10.0^y,  text("10^{$y}", 9, :left, RGB(0.75,0.75,0.75), :top)))
#     end

#     #annotate axes
#     # annotate!(f,(0.07*xlims(f)[2],1.8*ylims(f)[1],text(L"\textbf{Size}\rightarrow", 10, :left, RGB(0.25,0.25,0.25), :top)))
#     # annotate!(f,(xlims(f)[1],20*ylims(f)[1],text(L"\leftarrow\textbf{Conductance}", 20, :left, RGB(0.25,0.25,0.25), :top)))
  
#     return f
# end

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
        nodes = collect(union(sets[s]...))
        bin_val = (sum(snodes[nodes])/length(nodes))/ntrials 
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




function makefigs(gname; title="", nlines = 3, 
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
    
    ydim = (nlines+1)*35 + (isempty(title) ? 0 : 15) + (nlines == 2 ? 1 : 0)
    ylims = (10.0^(-nlines-1),2.0)
    #perform plotting 
    f = plot(margin=-20Measures.mm,  topmargin=isempty(title) ? 0Measures.mm : -nlines*0.5Measures.mm , size=(250,ydim), #background_color_inside=RGB(0.9,0.9,0.9),
                ylims=ylims, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
    makegridlines!(f, [1,2,3,4,5], -(1:nlines),tickfontsize=tickfontsize)
    
    #missed sets plots    
    tmpfs = _plot_missed_sets1!(deepcopy(f),gname,ncploc=ncploc,sloc=sloc,ntrials=50,model="seir",
                betas=betas,nbins=nbins,color=color)

    return tmpfs
end

function makefigs_1(gname; title="", nlines = 3, 
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
    
    ydim = (nlines+1)*35 + (isempty(title) ? 0 : 15) + (nlines == 2 ? 1 : 0)
    ylims = (10.0^(-nlines-1),2.0)
    @show ydim, ylims
    #perform plotting 
    f = plot(margin=-20Measures.mm,  topmargin=isempty(title) ? 0Measures.mm : -nlines*0.5Measures.mm , size=(250,ydim), #background_color_inside=RGB(0.9,0.9,0.9),
                ylims=ylims, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
    makegridlines!(f, [1,2,3,4,5], -(1:nlines),tickfontsize=tickfontsize)
    
    #missed sets plots    
    tmpfs = _plot_missed_sets1!(deepcopy(f),gname,ncploc=ncploc,sloc=sloc,ntrials=50,model="seir",
                betas=betas,nbins=nbins,color=color)

    return tmpfs
end

#generate individual plots 
gs = ["mexico","dblp","anon","study-11-2022-50","geometric"]
gs = vcat(getgnames.(gs,"input/graphs/")...)

fig_params = Dict() #gname->(betas,qpercents)
fig_params["modmexico-city.smat"] = collect(2e-2:2e-2:1e-1)
fig_params["dblp-cc.smat"] = collect(2e-2:2e-2:1e-1)
fig_params["anony-interactions-onemonthA-cc.smat"] = collect(2e-2:2e-2:1e-1)
fig_params["study-11-2022-50.smat"] = vcat(4e-3:2e-3:1e-2,2e-2)
fig_params["geometric-100000-2d-1-20.smat"] = vcat(7e-3:1e-3:1e-2,2e-2)


c = cgrad(cgrad(:magma)[1:240])

#do each plot manually
gname = gs[1]
nlines=4
nbins = 55
tmpfigs = makefigs(gname,nlines=nlines,betas=fig_params[gname],color=c,nbins=nbins)
fs = [tmpfigs[(beta,qpercent)] for (beta,qpercent) in sort(collect(keys(tmpfigs)))]
plot!.(fs,title="",ylims=(3e-5,1),xlims=(1,2e5),size=(600,400),clims=(0,1))

f = make_acl_fig(gname,nlines=nlines,nbins=nbins)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))

newf = plot(vcat([f],fs)...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)
plot!(newf,dpi=500)
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets-$(gname[1:end-5]).png")


gname = gs[2]
nlines = 3
nbins = 60
tmpfigs = makefigs(gname,nlines=nlines,betas=fig_params[gname],color=c,nbins=nbins)
fs = [tmpfigs[(beta,qpercent)] for (beta,qpercent) in sort(collect(keys(tmpfigs)))]
plot!.(fs,title="",ylims=(3e-4,1),xlims=(1,2e5),size=(600,400),clims=(0,1))

f = make_acl_fig(gname,nlines=nlines,nbins=nbins)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))

newf = plot(vcat([f],fs)...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)

plot!(newf,dpi=500)
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets-$(gname[1:end-5]).png")


gname = gs[3]
nlines = 3
nbins = 60
tmpfigs = makefigs(gname,nlines=nlines,betas=fig_params[gname],color=c,nbins=nbins)
fs = [tmpfigs[(beta,qpercent)] for (beta,qpercent) in sort(collect(keys(tmpfigs)))]
plot!.(fs,title="",ylims=(3e-4,1),xlims=(1,2e5),size=(600,400),clims=(0,1))

f = make_acl_fig(gname,nlines=nlines,nbins=nbins)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))

newf = plot(vcat([f],fs)...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)
plot!(newf,dpi=500)
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets-$(gname[1:end-5]).png")
            

gname = gs[4]
nlines = 3
nbins = 60
tmpfigs = makefigs(gname,nlines=nlines,betas=fig_params[gname],color=c,nbins=nbins)
fs = [tmpfigs[(beta,qpercent)] for (beta,qpercent) in sort(collect(keys(tmpfigs)))]
plot!.(fs,title="",ylims=(3e-4,1),xlims=(1,2e5),size=(600,400),clims=(0,1))

f = make_acl_fig(gname,nlines=nlines,nbins=nbins)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))

newf = plot(vcat([f],fs)...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)
plot!(newf,dpi=500)
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets-$(gname[1:end-5]).png")
                        


#getting high quality image of colorbar for mexico city
# gname = "modmexico-city.smat"
# tmpfigs = makefigs_1(gname,nlines=4,betas=fig_params[gname],color=c,nbins=50)
# fs = [tmpfigs[(beta,qpercent)] for (beta,qpercent) in sort(collect(keys(tmpfigs)))]
# plot!.(fs,title="",ylims=(1.e-5,1.15),xlims=(0.8,2e5),size=(600,400),clims=(0,1))

# f = deepcopy(fs[2])
# plot!(f,colorbar=true,clims=(0,1),
#         bottom_margin = -5*Measures.mm,
#         left_margin=-14Measures.mm,
#         dpi=1000)
# savefig(f,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/mexico-city-colorbar.png")



function makefigs1(gname; title="", nlines = 3, 
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
    
    ydim = (nlines+1)*35 + (isempty(title) ? 0 : 15) + (nlines == 2 ? 1 : 0)
    ylims = (10.0^(-nlines-1),2.0)
    @show ydim, ylims
    #perform plotting 
    f = plot(margin=-20Measures.mm,  topmargin=isempty(title) ? 0Measures.mm : -nlines*0.5Measures.mm , size=(250,ydim), #background_color_inside=RGB(0.9,0.9,0.9),
                ylims=ylims, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
    makegridlines!(f, [1,2,3,4,5], -(1:nlines),tickfontsize=tickfontsize)
    
    #missed sets plots    
    tmpfs = _plot_missed_sets1!(deepcopy(f),gname,ncploc=ncploc,sloc=sloc,ntrials=50,model="seir",
                betas=betas,nbins=nbins,color=color)

    return tmpfs
end

function make_missed_sets_plot(gname::String,betas::Vector{Float64};
                nlines::Int=4,nbins::Int=60,color=:inferno)

    figs = makefigs1(gname,nlines=nlines,betas=betas,color=color,nbins=nbins)

    fs = [figs[(beta,qpercent)] for (beta,qpercent) in sort(collect(keys(figs)))]
    plot!.(fs,title="",xlims=(1,2e5),size=(600,400),clims=(0,1))

    return fs
end

function preview_plots(fs)
    newf = plot(fs...,layout = (lastindex(fs),1),
        size=(650,500*lastindex(fs)),
        top_margin=1Measures.mm,
        bottom_margin=0Measures.mm,
        left_margin=0Measures.mm,
        right_margin=0Measures.mm)

    return newf
end


#doing other graphs from the paper 
hs = ["commutes-all", "filtered", #row1
    "cn-moduillinois", "cn-Penn", "cn-modWiscon", #row2...
    "enron", #row 3
    "cit-HepPh", "slashdot", #"flickr", #skip flickr for now 
    "geometric","cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.9-0.0-100-5.smat"]
    #  "rewired-10000.0-modmexico","rewired-10000.0-cn-moduillinois","er-10000.0-dblp-cc"]
gnames = map(h->h=="" ? "" : getgnames(h,"pipeline/graphs/")[1], hs)
c = cgrad(cgrad(:magma)[1:240])

#### Manually doing graphs.. need to inspect them and adjust 

master_figs = []
#commutes
gname = "commutes-all.smat"
fs = make_missed_sets_plot(gname,get_betas(gname),nlines=2,color=c)
inds = 1:1:5

f = make_acl_fig(gname,nlines=2,nbins=60)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))
newf = plot(vcat([f],fs[inds])...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)
push!(master_figs,deepcopy(newf))
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets-$(gname[1:end-5])")
get_betas(gname)[inds]


#filtered us flows 
gname = "covidflows-2020_08_31-filtered-20.smat"
fs = make_missed_sets_plot(gname,get_betas(gname),nlines=3,color=c)
# preview_plots(fs)
inds = 13:2:21
f = make_acl_fig(gname,nlines=3,nbins=60)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))
newf = plot(vcat([f],fs[inds])...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)

push!(master_figs,deepcopy(newf))
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets-$(gname[1:end-5])")
get_betas(gname)[inds]


#cn-uill
gname = "cn-modUIllinois20.smat"
fs = make_missed_sets_plot(gname,get_betas(gname),nlines=2,color=c)
# preview_plots(fs)
inds = 13:3:25
f = make_acl_fig(gname,nlines=2,nbins=60)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))
newf = plot(vcat([f],fs[inds])...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)
push!(master_figs,deepcopy(newf))
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets-$(gname[1:end-5])")
get_betas(gname)[inds]


#cn-penn
gname = "cn-Penn94.smat"
fs = make_missed_sets_plot(gname,get_betas(gname),nlines=2,color=c)
preview_plots(fs)
inds = 13:3:25
f = make_acl_fig(gname,nlines=2,nbins=60)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))
newf = plot(vcat([f],fs[inds])...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)

push!(master_figs,deepcopy(newf))
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets-$(gname[1:end-5])")
get_betas(gname)[inds]


#cn-wisc
gname = "cn-modWisconsin87.smat"
fs = make_missed_sets_plot(gname,get_betas(gname),nlines=2,color=c)
# preview_plots(fs)
inds= 13:3:25
f = make_acl_fig(gname,nlines=2,nbins=60)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))
newf = plot(vcat([f],fs[inds])...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)

push!(master_figs,deepcopy(newf))
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets-$(gname[1:end-5])")
get_betas(gname)[inds]

#enron
gname = "email-Enron.smat"
fs = make_missed_sets_plot(gname,get_betas(gname),nlines=3,color=c)
# preview_plots(fs)
inds = 6:2:14
f = make_acl_fig(gname,nlines=3,nbins=60)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))
newf = plot(vcat([f],fs[inds])...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)

push!(master_figs,deepcopy(newf))
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets-$(gname[1:end-5])")
get_betas(gname)[inds]

#citation
gname = "cit-HepPh.smat"
fs = make_missed_sets_plot(gname,get_betas(gname),nlines=2,color=c)
# preview_plots(fs)
inds = 4:2:12
f = make_acl_fig(gname,nlines=2,nbins=60)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))
newf = plot(vcat([f],fs[inds])...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)

push!(master_figs,deepcopy(newf))
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets-$(gname[1:end-5])")
get_betas(gname)[inds]

#slashdot
gname = "Slashdot0811.smat"
fs = make_missed_sets_plot(gname,get_betas(gname),nlines=2,color=c)
# preview_plots(fs)
inds = 3:2:11
f = make_acl_fig(gname,nlines=2,nbins=60)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))
newf = plot(vcat([f],fs[inds])...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)

push!(master_figs,deepcopy(newf))
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets-$(gname[1:end-5])")
get_betas(gname)[inds]

#local geometric 
gname = "geometric-100000-2d-1-20.smat"
fs = make_missed_sets_plot(gname,get_betas(gname),nlines=3,color=c)
# preview_plots(fs)
inds = 8:1:12
f = make_acl_fig(gname,nlines=3,nbins=60)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))
newf = plot(vcat([f],fs[inds])...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)

push!(master_figs,deepcopy(newf))
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets-$(gname[1:end-5])")
get_betas(gname)[inds]

#random walk communities
gname = "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.9-0.0-100-5.smat"
fs = make_missed_sets_plot(gname,get_betas(gname),nlines=4,color=c)
# preview_plots(fs)

inds = 9:2:17
f = make_acl_fig(gname,nlines=4,nbins=60)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))
newf = plot(vcat([f],fs[inds])...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)

push!(master_figs,deepcopy(newf))
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets-$(gname[1:end-5]).png")
get_betas(gname)[inds]

#redoing some of the old ones for a potentially complete figure
#mexico
gname = "modmexico-city.smat"
fs = make_missed_sets_plot(gname,get_betas(gname),nlines=4,color=c)
# preview_plots(fs)

inds = 11:2:19
f = make_acl_fig(gname,nlines=4,nbins=60)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))
newf = plot(vcat([f],fs[inds])...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)

push!(master_figs,deepcopy(newf))
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets1-$(gname[1:end-5]).png")
get_betas(gname)[inds]

#dblp
gname = "dblp-cc.smat"
fs = make_missed_sets_plot(gname,get_betas(gname),nlines=3,color=c)
# preview_plots(fs)

inds = 11:2:19
f = make_acl_fig(gname,nlines=3,nbins=60)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))
newf = plot(vcat([f],fs[inds])...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)

push!(master_figs,deepcopy(newf))
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets1-$(gname[1:end-5]).png")
get_betas(gname)[inds]

#facebook interactions
gname = "anony-interactions-onemonthA-cc.smat"
inds = 11:2:19

fs = make_missed_sets_plot(gname,get_betas(gname),nlines=3,color=c)

f = make_acl_fig(gname,nlines=3,nbins=60)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))
newf = plot(vcat([f],fs[inds])...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)

push!(master_figs,deepcopy(newf))
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets1-$(gname[1:end-5]).png")
get_betas(gname)[inds]

#geometric communities
gname = "study-11-2022-50.smat"
fs = make_missed_sets_plot(gname,get_betas(gname),nlines=3,color=c)
# preview_plots(fs)

inds = 7:1:11
f = make_acl_fig(gname,nlines=3,nbins=60)
plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))
newf = plot(vcat([f],fs[inds])...,layout=(1,6),size=(600*6,400),
            left_margin=-15Measures.mm,
            right_margin=-3Measures.mm,
            bottom_margin=-3Measures.mm,
            top_margin=0Measures.mm)

push!(master_figs,deepcopy(newf))
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/missed-sets1-$(gname[1:end-5]).png")
get_betas(gname)[inds]


master_figs

Plots.get_size(master_figs[2])
plot(master_figs...,layout=(14,1),size = (3700,450*14),
        bottom_margin=2Measures.mm,
        top_margin=2Measures.mm,
        right_margin=1Measures.mm)

preview_plots(master_figs)
## trying to dig deeper into issues with flickr..


ylims(master_figs[1][2])

#generating a figure that explains what we are seeing here
#(acl ncp) - (epi ncp) = (missed sets ncp)

#ncp info for both processes
# epidemic ncp
# 

######################################################################
################ Panels for Explaining Missed Sets ###################
######################################################################

#testing
# xbounds,ybounds = (2.5,120000),(5e-5,1)
# f = plot(margin=-20Measures.mm,  topmargin=isempty("") ? 0Measures.mm : -nlines*0.5Measures.mm , size=(900,600), 
#               ylims=ybounds, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
# makegridlines!(f, [1,2,3,4,5], -(1:4),tickfontsize=15)

gname = "modmexico-city.smat"

using ColorSchemes
### acl ncp figure
#load data 


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


function load_acl_ncp(gname;datadir="pipeline/data/",normalizex::Bool=true)
    base_graph = canonical_graph_name(gname)
    
    ncp,headerinfo = readdlm(joinpath(datadir,base_graph[1:end-5],"ncpdata","ncpinfo-$(gname[1:end-5]).txt"),',',header=true)
    ncp = DataFrame(ncp,vec(headerinfo))
    
    # x = map(x->min(ncp.size[x],n-ncp.size[x]),1:size(ncp,1))
    if normalizex 
        n = graphSize(gname)
        x = map(x->min(ncp.size[x],n-ncp.size[x]),1:size(ncp,1))
    else 
        x = ncp.size
    end

    y = ncp.cond

    return x,y
end 

function acl_ncp_fig(gname::String;color=:cividis,nbins=80,nlines=4)
    x,y = load_acl_ncp(gname)
    xbounds,ybounds = (1.0,2e5),(1.52e-5,1.15)
    f = plot(margin=-20Measures.mm,  topmargin=isempty("") ? 0Measures.mm : -nlines*0.5Measures.mm , size=(600,400), 
                  ylims=ybounds, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
    makegridlines!(f, [1,2,3,4,5], -(1:nlines),tickfontsize=15)

    myhexbin2!(f,x,y,xlims=xbounds,ylims=ybounds,color=color,nbins=nbins,normalize_color=true)
    myhexbin2!(f,x,y,xlims=xbounds,ylims=ybounds,color=color,nbins=nbins,normalize_color=true)
    return f
end

# acl_color = cgrad(cgrad(:viridis)[1:250])
# f = acl_ncp_fig(gname,color=acl_color,nbins=55,nlines=4)
# plot!(f,colorbar=true,
#         bottom_margin=-5Measures.mm,
#         right_margin=-1Measures.mm,
#         dpi=1000)
# savefig(f,"/p/mnt/scratch/network-epi/code/paper-figs/missed-sets-figs/acl-ncp-$(gname[1:end-5]).png")



# function acl_ncp_fig(gname::String;color=:heat,nbins=80)
#     x,y = load_acl_ncp(gname)
#     xbounds,ybounds = (1.0,2e5),(3e-5,1.0)
#     f = plot(margin=-20Measures.mm,  topmargin=isempty("") ? 0Measures.mm : -nlines*0.5Measures.mm , size=(600,400), 
#                   ylims=ybounds, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
#     makegridlines!(f, [1,2,3,4,5], -(1:4),tickfontsize=15)

#     myhexbin2!(f,x,y,
#             # xlims=xbounds,ylims=ybounds,
#             color=color,nbins=nbins,normalize_color=true)
#     myhexbin2!(f,x,y,
#             # xlims=xbounds,ylims=ybounds,
#             color=color,nbins=nbins,normalize_color=true)
#     return f
# end


# ydim = (nlines+1)*35 + (isempty(title) ? 0 : 15) + (nlines == 2 ? 1 : 0)
# ylims = (10.0^(-nlines-0.5),2.0)
# @show ydim, ylims
# #perform plotting 
# f = plot(margin=-20Measures.mm,  topmargin=isempty(title) ? 0Measures.mm : -nlines*0.5Measures.mm , size=(250,ydim), #background_color_inside=RGB(0.9,0.9,0.9),
#             ylims=ylims, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
# makegridlines!(f, [1,2,3,4,5], -(1:nlines),tickfontsize=tickfontsize)


function make_acl_fig(gname; title="", nlines = 3, 
                    # ybounds = (1e-3,1),
                    xbounds = (0.8,2e5),
                    tickfontsize=15,
                    nbins=60,
                    color=cgrad(cgrad(:viridis)[1:250]))
    
    ydim = (nlines+1)*35 + (isempty(title) ? 0 : 15) + (nlines == 2 ? 1 : 0)
    ylims = (10.0^(-nlines-1),2.0)
    @show ydim, ylims
    #perform plotting 
    f = plot(margin=-20Measures.mm,  topmargin=isempty(title) ? 0Measures.mm : -nlines*0.5Measures.mm , size=(250,ydim), #background_color_inside=RGB(0.9,0.9,0.9),
                ylims=ylims, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
                makegridlines!(f, [1,2,3,4,5], -(1:nlines),tickfontsize=tickfontsize)

    x,y = load_acl_ncp(gname)
    myhexbin2!(f,x,y,
            # xlims=xbounds,ylims=ybounds,
            color=color,nbins=nbins,normalize_color=true)
    myhexbin2!(f,x,y,
            # xlims=xbounds,ylims=ybounds,
            color=color,nbins=nbins,normalize_color=true)
    return f
end

          



#gname, nlines, nbins, betas,
function make_newplot(gname,betas,nlines::Int=3,nbins=60)
    c = cgrad(cgrad(:magma)[1:240])
    tmpfigs = makefigs(gname,nlines=nlines,betas=betas,color=c,nbins=nbins,xbounds=(0.8,2e5))
    fs = [tmpfigs[(beta,qpercent)] for (beta,qpercent) in sort(collect(keys(tmpfigs)))]
    plot!.(fs,title="",ylims=(3e-5,1),xlims=(0.8,2e5),size=(600,400),clims=(0,1))

    f = make_acl_fig(gname,nlines=nlines,nbins=nbins)
    plot!(f,xlims=xlims(fs[1]),ylims=ylims(fs[1]))

    # newf = plot(vcat([f],fs)...,layout=(1,6),size=(600*6,400),
    #         left_margin=-15Measures.mm,
    #         right_margin=-3Measures.mm,
    #         bottom_margin=-3Measures.mm,
    #         top_margin=0Measures.mm)
    return vcat([f],fs)
end


#parameter list 
#gname,betas,nlines
hs = ["commutes-all", "mexico","filtered", #row1
    "cn-moduillinois", "cn-Penn", "cn-modWiscon", #row2...
    "dblp", "enron", "anon", #row 3
    "cit-HepPh", "slashdot", #"flickr", #skip flickr for now 
    "geometric","study-11-2022-50","cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.9-0.0-100-5.smat"]
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
    vcat(4e-3:2e-3:1e-2,2e-2), #geo-comm 
    [0.009, 0.02, 0.04, 0.06, 0.08], #rwc 
]

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
    3, #geo-comm  
    4, #rwc 
]


fvec = make_newplot(gnames[1],beta_vec[1],nlines_vec[1])
plot!.(fvec,ylims=(1e-3,1.1),size=(600,400))
fvec[1]
fvec[2]
xlims(fvec[1])
xlims(fvec[2])