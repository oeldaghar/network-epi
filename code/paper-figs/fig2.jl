parent_path = "/p/mnt/scratch/network-epi/"
include(joinpath(parent_path,"code/graph-io.jl")) 

## plotting fcns
include(joinpath(parent_path,"code/ncp/hexbins1.jl"))
# include("ncpplots.jl")
include(joinpath(parent_path,"code/ncp/ncpplots1.jl"))

ENV["GKSwstype"] = "100"
gpath = "pipeline/graphs/"


include(joinpath(parent_path,"code/fast-diffusion.jl"))
include(joinpath(parent_path,"code/ncp/parallel-epidemic-ncp.jl"))


using Measures
using Plots

gpath = "pipeline/graphs/"
gs = readdir("input/graphs/")
filter!(x->endswith(x,".smat"),gs)
#ignoring these for now 
filter!(c->!occursin("flickr",lowercase(c)),gs)
filter!(c->!occursin("livejournal",lowercase(c)),gs)

total_trials = 60000
#horizontal then vertical 
hs = ["dblp","enron","",
    "anon", "", "",
    "cn-moduillinois", "cn-Penn", "cn-modWiscon",
    "commutes-all","mexico", "",
    "geometric","study-11-45","cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.9-0.0-100-5.smat",
     "rewired-10000.0-modmexico","rewired-10000.0-cn-moduillinois","er-10000.0-dblp-cc"]


gnames = []
for h in hs
    if h==""
        push!(gnames,h)
    else
        push!(gnames,getgnames(h,gpath)[1])
    end
end

function makegridlines!(f,xs,ys)
    #jerry rig this to work for now. 
    #cant seem to find in the docs...

    #vertical lines
    for x in xs
        plot!(f,[x;x],[ylims(f)...],c=:grey,linewidth=1,leg=false,alpha=0.8)
    end

    #horizontal lines
    for y in ys
        plot!(f,[xlims(f)...],[y;y],c=:grey,linewidth=1,leg=false,alpha=0.8)
        # hline!(f,[1e-1,1e-2,1e-3],label="", color=:grey, linewidth=1)
    end
    return f
end

xlimsoffset = 2
xlabeloffset = 0.5
glowxoffset = .15
function myxlabel!(f,xval, x, y, str, args...; kwargs...) 
    #miny = minimum(xy->xy[2], filter(xy->abs(log(xy[1])-log(xval)) <= 1, collect(zip(x,y))))
    miny,maxy = extrema(xy->xy[2], filter(xy->abs(log(xy[1])-log(xval)) <= 1, collect(zip(x,y))))
    plot!(f,[xval,xval],[10*miny, miny*xlabeloffset], linecolor=:grey, label="")
    # very hacky glow 
    for i in 1:100
        annotate!(f,[(xval*(1 + glowxoffset*(2*rand()-1)), 
        miny*xlabeloffset*(1 + glowxoffset*(2*rand()-1)), 
        text(str, :white; pointsize=11, 
        kwargs..., valign=:top, ))])
    end 
    annotate!(f,[(xval, miny*xlabeloffset, text(str, :black; pointsize=12, 
        kwargs..., valign=:top, ))]) 
end 

#tick labels 
offsetscale = 0.3 # pick this to get the text right on the line... 
tlabel = (f,y,str) -> annotate!(f,[(xlims()[2], y*offsetscale, text(str, :black; pointsize=12, 
  halign=:right, valign=:bottom, ))])


#main fcn for generating individual plots 
function make_fig(gnames,gpath="pipeline/graphs/")
    #supress all axis not on the edge 
    #make raw plots and save in the same shape as gnames 
    xbounds = (0.8,200000)
    ybounds = (5e-4,1.5)
    figs = []
    for gname in gnames
        #load in data 
        f = plot(leftmargin=-10mm, topmargin=-0mm, bottommargin=-10mm, rightmargin=-10mm, size=(250,175),
            ylims=(1e-4,1), framestyle=:none, xlims=(0.5,xlimsoffset*1e5))

        if gname!=""
            #load data 
            A = loadGraph(joinpath(gpath,gname))
            if startswith(gname,"er-10000.0") 
                dst = joinpath("pipeline/data/$(gname[12:end-5])/ncpdata/")    
            elseif startswith(gname,"rewired-10000.0")
                dst = joinpath("pipeline/data/$(gname[17:end-5])/ncpdata/")
            else
                dst = joinpath("pipeline/data/$(gname[1:end-5])/ncpdata/")
            end
            
            ncp,headerinfo = readdlm(joinpath(dst,"ncpinfo-epidemic-seir-60000-$(gname[1:end-5]).txt"),',',header=true)
            ncp = DataFrame(ncp,vec(headerinfo))
            
            x = map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1))
            y = ncp.cond
            
            #perform plotting 
            f = plot(leftmargin=-10mm, topmargin=-0mm, bottommargin=-10mm, rightmargin=-10mm, size=(500,375),
                        ylims=(1e-4,1), framestyle=:none, xlims=(0.5,xlimsoffset*1e5))
            makegridlines!(f,[],[1e-1,1e-2,1e-3])
            #add labels for mexicocity
            if gname=="modmexico-city.smat"
                #annotated xticks 
                myxlabel!(f,1e1, x, y, "10^{1}", halign=:center)
                myxlabel!(f,1e2, x, y, "10^{2}", halign=:center)
                myxlabel!(f,1e3, x, y, "10^{3}", halign=:center)
                myxlabel!(f,1e4, x, y, "10^{4}", halign=:center)
            end
            if gname=="commutes-all.smat"
                #annotated yticks 
                tlabel(f,1e-1, "10^{-1}")
                tlabel(f,1e-2, "10^{-2}")
                tlabel(f,1e-3, "10^{-3}")
            end
            #place ncp data on plot 
            f = myncpplot2!(f,x,y,
                xlims=xbounds,ylims=ybounds,nbins=(120,100))#
                # plotmin=false)
                # plotmedian=false)
            #title 
            if occursin("cl-lfr",gname)
                title!(f,"cl-lfr-rw-$(gname[72:75])")
            elseif occursin("anony-interactions",gname)
                title!(f,"anony-interactions")
            elseif occursin(Regex("rewired-[\\d\\.\\d]*"),gname)
                title!(f,"rewired-$(gname[17:end-5])")
            elseif occursin(Regex("er-[\\d\\.\\d]*"),gname)
                title!(f,"er-$(gname[12:end-5])")
            else
                title!(f,"$(gname[1:end-5])")
            end
        end

        push!(figs,f)
    end
    return figs
end


# Plots.plot(rand(10,12),layout=(4,3),xaxis=true)
figs = make_fig(gnames)

fs = deepcopy(figs)

#combine plots into single combined figure 
newf = Plots.plot(fs...,layout=(6,3),titlefont=font(12),
    size=(750,1200),margin=-2mm)#,thickness_scaling=0.8,tickfontsize=font(5))
plot!(newf,dpi=500)
Plots.savefig(newf,"code/paper-figs/fig2/epidemic-ncps-figure-2.png")






xbounds = (0.8,200000)
ybounds = (5e-4,1.5)
figs = []
gname = gnames[11]

#load in data 
f = plot(framestyle=:none,xlims=xbounds,ylims=ybounds,
        margins=-5mm,leg=false)


A = loadGraph(joinpath(gpath,gname))
if startswith(gname,"er-10000.0") 
    dst = joinpath("pipeline/data/$(gname[12:end-5])/ncpdata/")    
elseif startswith(gname,"rewired-10000.0")
    dst = joinpath("pipeline/data/$(gname[17:end-5])/ncpdata/")
else
    dst = joinpath("pipeline/data/$(gname[1:end-5])/ncpdata/")
end

ncp,headerinfo = readdlm(joinpath(dst,"ncpinfo-epidemic-seir-60000-$(gname[1:end-5]).txt"),',',header=true)
ncp = DataFrame(ncp,vec(headerinfo))
# f = myncpplot1(map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1)),ncp.cond,
#     plotmin=false,plotmedian=false,xlims=(1,200000),ylims=(1e-3,1.1))

#make gridlines for plots where there is data 
# f = makegridlines!(f,[1;10;100;1000;10000;100000],[1;1e-1;1e-2;1e-3])
f = makegridlines!(f,[],[1e-1;1e-2;1e-3])
#add ncpdata to plot 
x = map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1))
y = ncp.cond
if gname=="modmexico-city.smat"
    myxlabel!(f,1e1, x, y, "10^{1}", halign=:center)
    myxlabel!(f,1e2, x, y, "10^{2}", halign=:center)
    myxlabel!(f,1e3, x, y, "10^{3}", halign=:center)
    myxlabel!(f,1e4, x, y, "10^{4}", halign=:center)
end

f = myncpplot2!(f,x,y,
    xlims=xbounds,ylims=ybounds,nbins=(120,100))
#title 
if occursin("cl-lfr",gname)
    title!(f,"cl-lfr-rw-$(gname[72:75])")
elseif occursin("anony-interactions",gname)
    title!(f,"anony-interactions")
elseif occursin(Regex("rewired-[\\d\\.\\d]*"),gname)
    title!(f,"rewired-$(gname[17:end-5])")
elseif occursin(Regex("er-[\\d\\.\\d]*"),gname)
    title!(f,"er-$(gname[12:end-5])")
else
    title!(f,"$(gname[1:end-5])")
end
plot!(f,titlefontsize=5,margins=-5mm)

f = plot!(f,xlims=xbounds,ylims=ybounds,margins=-20mm)



#############################################################
ncp,headerinfo = readdlm(joinpath(dst,"ncpinfo-epidemic-seir-60000-$(gname[1:end-5]).txt"),',',header=true)
ncp = DataFrame(ncp,vec(headerinfo))


f = myncpplot2!(f,x,y,
    xlims=xbounds,ylims=ybounds,nbins=(120,100))
#title 
if occursin("cl-lfr",gname)
    title!(f,"cl-lfr-rw-$(gname[72:75])")
elseif occursin("anony-interactions",gname)
    title!(f,"anony-interactions")
elseif occursin(Regex("rewired-[\\d\\.\\d]*"),gname)
    title!(f,"rewired-$(gname[17:end-5])")
elseif occursin(Regex("er-[\\d\\.\\d]*"),gname)
    title!(f,"er-$(gname[12:end-5])")
else
    title!(f,"$(gname[1:end-5])")
end
plot!(f,margins=-5mm)
tlabel(f,1e-1, "10^{-1}")
tlabel(f,1e-2, "10^{-2}")
tlabel(f,1e-3, "10^{-3}")  
f


