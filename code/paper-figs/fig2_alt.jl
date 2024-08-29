#allows for execution from command line as well as an ide so files can be run modularlly
mainDir = joinpath("/",split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))]...)

include(joinpath(mainDir,"code/graph-io.jl")) 
include(joinpath(mainDir,"code/ncp/ncpplots1.jl")) 

using DelimitedFiles
using DataFrames
using Measures
using Pkg
using PerceptualColourMaps
using LaTeXStrings

##
gr()
# ENV["GRDIR"] = ""
# Pkg.build("GR")

function graphSize(gname)
  A = loadGraph(joinpath("pipeline/graphs/",gname))
  return size(A,1)
end

##
function load_ncp(gname,datadir="pipeline/data/",normalizex::Bool=true)
  base_graph = canonical_graph_name(gname)
  
  ncp,headerinfo = readdlm(joinpath(datadir,base_graph[1:end-5],"ncpdata","ncpinfo-epidemic-subsampling-4-seir-50000-$(gname[1:end-5]).txt"),',',header=true)
  ncp = DataFrame(ncp,vec(headerinfo))
  
  # x = map(x->min(ncp.size[x],n-ncp.size[x]),1:size(ncp,1))
  if normalizex 
    n = graphSize(gname)
    x = map(x->min(ncp.size[x],n-ncp.size[x]),1:size(ncp,1))
  else 
    x = ncp.size
  end
  # x = ncp.size
  y = ncp.cond

  return x,y
end 

function load_ncp_cut(gname,datadir="pipeline/data/",normalizex::Bool=true)
  base_graph = canonical_graph_name(gname)
  
  ncp,headerinfo = readdlm(joinpath(datadir,base_graph[1:end-5],"ncpdata","ncpinfo-epidemic-subsampling-4-seir-50000-$(gname[1:end-5]).txt"),',',header=true)
  ncp = DataFrame(ncp,vec(headerinfo))
  
  # x = map(x->min(ncp.size[x],n-ncp.size[x]),1:size(ncp,1))
  if normalizex 
    n = graphSize(gname)
    x = map(x->min(ncp.size[x],n-ncp.size[x]),1:size(ncp,1))
  else 
    x = ncp.size
  end
  # x = ncp.size
  y = ncp.cut

  return x,y
end

#testing  
# gname = getgnames("sbm","input/graphs/")[1]
# x,y = load_ncp(gname)

# f = plot()
# myhexbin_conditional!(f,x,y;
#         nbins=(80,80),
#         color=cgrad(cmap("CBL2")[75:230])
#         )
# plot!(f,colorbar=true,tickfontsize=10,dpi=500,clims=(0,1))



# savefig(f, "code/paper-figs/fig2/large-colorbar.png")
# myncpplot_conditional!(f,x,y,plotmin=false)
# extrema(x)

# A = loadGraph("cit-Patents.smat","input/graphs/")
# myncpplot_conditional!(f,x,y,xlims=extrema(x),ylims=extrema(y))


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



function makefig(gname; title="", nlines = 3, ybounds = (1e-3,1), xbounds = (0.8,2e5), nbins=(80,60),
                  color=:inferno,clims=(0,1),tickfontsize=10)
  println(gname)
  ydim = (nlines+1)*35 + (isempty(title) ? 0 : 15) + (nlines == 2 ? 1 : 0)
  ylims = (10.0^(-nlines-1),2.0)
  @show ydim, ylims
  #perform plotting 
  f = plot(margin=-20Measures.mm,  topmargin=isempty(title) ? 0Measures.mm : -nlines*0.5Measures.mm , size=(250,ydim), #background_color_inside=RGB(0.9,0.9,0.9),
              ylims=ylims, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
  
  
  if gname != "" 
    #makegridlines!(f, [1e1,1e2,1e3,1e4,1e5], map(i->10.0^(-i), 1:nlines))
    makegridlines!(f, [1,2,3,4,5], -(1:nlines),tickfontsize=tickfontsize)
    x,y = load_ncp(gname)

    #place ncp data on plot 
    myhexbin_conditional!(f,x,y;
        xlims=xbounds,ylims=ybounds,nbins=nbins,
        color=color)
    #plotting a second time 
    myhexbin_conditional!(f,x,y;
        xlims=xbounds,ylims=ybounds,nbins=nbins,
        color=color)
    # myncpplot_conditional!(f,x,y;
    #     xlims=xbounds,ylims=ybounds,nbins,#
    #     plotmin=false, markersize=0, linecolor=:orange, linewidth=2,
    #     plotmedian=false,color=color)

    plot!(f;ylims=ylims,clims=clims)
    
    title!(f,title, 
              titlelocation=:center, titlevalign=:top,
              titlefontsize=15, titlefont="Helvetica")
  else
    title!(f," ", 
              titlelocation=:left, titlevalign=:top,
              titlefontsize=15, titlefont="Helvetica")
  end
  return f
end

function makefig_cut(gname; title="", nlines = 6, ybounds = (1,1e7), xbounds = (0.8,2e5), nbins=(80,60),
                  color=:inferno,clims=(0,1),tickfontsize=10)
  println(gname)
  ydim = (nlines+1)*35 + (isempty(title) ? 0 : 15) + (nlines == 2 ? 1 : 0)
  ylims = (0.95,2*ybounds[2])
  @show ydim, ylims
  #perform plotting 
  f = plot(margin=-20Measures.mm,  topmargin=isempty(title) ? 0Measures.mm : -nlines*0.5Measures.mm , size=(250,ydim), #background_color_inside=RGB(0.9,0.9,0.9),
              ylims=ylims, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
  
  
  if gname != "" 
    #makegridlines!(f, [1e1,1e2,1e3,1e4,1e5], map(i->10.0^(-i), 1:nlines))
    makegridlines!(f, [1,2,3,4,5], (1:nlines),tickfontsize=tickfontsize)
    x,y = load_ncp_cut(gname)

    #place ncp data on plot 
    myhexbin_conditional!(f,x,y;
        xlims=xbounds,ylims=ybounds,nbins=nbins,
        color=color)
    #plotting a second time 
    myhexbin_conditional!(f,x,y;
        xlims=xbounds,ylims=ybounds,nbins=nbins,
        color=color)
    # myncpplot_conditional!(f,x,y;
    #     xlims=xbounds,ylims=ybounds,nbins,#
    #     plotmin=false, markersize=0, linecolor=:orange, linewidth=2,
    #     plotmedian=false,color=color)

    plot!(f;ylims=ylims,clims=clims)
    
    title!(f,title, 
              titlelocation=:center, titlevalign=:top,
              titlefontsize=15, titlefont="Helvetica")
  else
    title!(f," ", 
              titlelocation=:left, titlevalign=:top,
              titlefontsize=15, titlefont="Helvetica")
  end
  return f
end

function _makefig1(gname; title="", nlines = 3, ybounds = (1e-3,1), xbounds = (0.8,2e5), nbins=(80,60),
                  color=:inferno,clims=(0,1),tickfontsize=10)
  println(gname)
  ydim = (nlines+1)*35 + (isempty(title) ? 0 : 15) + (nlines == 2 ? 1 : 0)
  ylims = (10.0^(-nlines-1),2.0)
  @show ydim, ylims
  #perform plotting 
  f = plot(margin=-20Measures.mm,  topmargin=isempty(title) ? 0Measures.mm : -nlines*0.5Measures.mm , size=(250,ydim), #background_color_inside=RGB(0.9,0.9,0.9),
              ylims=ylims, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
  
  
  
  if gname != "" 
    #makegridlines!(f, [1e1,1e2,1e3,1e4,1e5], map(i->10.0^(-i), 1:nlines))
    makegridlines!(f, [1,2,3,4,5], -(1:nlines),tickfontsize=tickfontsize)
    x,y = load_ncp(gname)

    #split x and y using set size 
    nnodes = graphSize(gname)
    midpt = round(nnodes/2)
    inds = findall(x.>midpt)
    other_inds = setdiff(1:lastindex(x),inds)

    x_missed = x[inds]
    y_missed = y[inds]

    x,y = x[other_inds],y[other_inds]
    
    #place ncp data on plot 
    myhexbin_conditional!(f,x,y;
        xlims=xbounds,ylims=ybounds,nbins=nbins,
        color=color)
    #plotting a second time 
    myhexbin_conditional!(f,x,y;
        xlims=xbounds,ylims=ybounds,nbins=nbins,
        color=color)

    #plotting missed sets using a different color 
    x_missed = map(a->min(a,nnodes-a),x_missed)
    myhexbin_conditional!(f,x_missed,y_missed;
        xlims=xbounds,ylims=ybounds,nbins=nbins,
        color=cgrad(cgrad(:inferno)[1:100]))
    
    plot!(f;ylims=ylims,clims=clims)
    
    title!(f,title, 
              titlelocation=:center, titlevalign=:top,
              titlefontsize=15, titlefont="Helvetica")
  else
    title!(f," ", 
              titlelocation=:left, titlevalign=:top,
              titlefontsize=15, titlefont="Helvetica")
  end
  #vertical line showing what happens once conductance flips 
  return f
end

# makefig("dblp-cc.smat"; title="Test", nlines=3,color=cgrad(cmap("L04")))

# makefig("flickr-links-sym.smat"; title="Flickr", nlines=3,color=cgrad(cmap("D04")))

# f = makefig("modmexico-city.smat"; title="cit-hep", nlines=3,color=cgrad(cmap("CBL2")[75:230]),clims=(0,1))
# plot!(f,size=(500,300))

#for ncp explanation figure
# getgnames("study","input/graphs/")[end]
f = makefig("study-20-draft-150.smat"; title="", nlines=2,color=cgrad(cmap("CBL2")[75:230]),
          xbounds=(0.8,550))
plot!(f,size=(400,250),dpi=1000,top_margin=-20*Measures.mm)
# Plots.savefig(f,"code/paper-figs/example-figs/epi-ncp-study-20-draft-150.png")


f = makefig("study-25-1.smat"; title="", nlines=2,color=cgrad(cmap("CBL2")[75:230]),
          xbounds=(0.8,40000))
plot!(f,size=(400,250),dpi=500,top_margin=-20*Measures.mm)
Plots.savefig(f,"code/paper-figs/spatial-figs/epi-ncp-study-25-1.png")

f = makefig("study-25-2.smat"; title="", nlines=2,color=cgrad(cmap("CBL2")[75:230]),
          xbounds=(0.8,40000))
plot!(f,size=(400,250),dpi=500,top_margin=-20*Measures.mm)
Plots.savefig(f,"code/paper-figs/spatial-figs/epi-ncp-study-25-2.png")

f = makefig("study-25-150.smat"; title="", nlines=2,color=cgrad(cmap("CBL2")[75:230]),
          xbounds=(0.8,40000))
plot!(f,size=(400,250),dpi=500,top_margin=-20*Measures.mm)
Plots.savefig(f,"code/paper-figs/spatial-figs/epi-ncp-study-25-150.png")


function cut_ncp_fig(gname::String,xbounds=(0.8,3e4),ybounds=(1e-3,1.2))
  x,y = load_ncp_cut(gname)
  f = myhexbin_conditional(x,y,nbins=(80,60),
    color=cgrad(cmap("CBL2")[75:230]),
    ylims=ybounds,
    xlims=xbounds)
  myhexbin_conditional!(f,x,y,nbins=(80,60),
    color=cgrad(cmap("CBL2")[75:230]),
    ylims=ybounds,
    xlims=xbounds)
  plot!(f,yticks=([1,10,100,1000,1e4,1e5,1e6]),
    xlabel="Size",ylabel="Size of Cut")
  return f
end

function ncp_fig(gname::String,xbounds=(0.8,3e4),ybounds=(1e-3,1.2);datadir="pipeline/data/")
  x,y = load_ncp(gname,datadir)
  f = myhexbin_conditional(x,y,nbins=(80,60),
    color=cgrad(cmap("CBL2")[75:230]),
    ylims=ybounds,
    xlims=xbounds)
  myhexbin_conditional!(f,x,y,nbins=(80,60),
    color=cgrad(cmap("CBL2")[75:230]),
    ylims=ybounds,
    xlims=xbounds)
  plot!(f,yticks=([1e-3,1e-2,1e-1,1]),
    xlabel="Size",ylabel="Conductance")
  return f
end

gname = getgnames("nested","input/graphs/")[1]
f = ncp_fig(gname,(0.8,5e4),(1e-4,1.2),datadir="/p/mnt/scratch/network-epi/pipeline/data/")
plot!(f,title="$(gname[1:end-5])\nnnodes: $(50*2^10), groupsize: $(split(gname,"-")[3]), layers: $(split(gname,"-")[4]),\navgd: $(split(gname,"-")[5]), decay_rate: $(split(gname,"-")[6])",
topmargin=6Measures.mm)

gnames = getgnames("study-25","input/graphs/")
gname = gnames[1]
ncp_fig(gname,(0.8,6e4))
plot!(title=gname,ylabel="Conductance")

using ProgressMeter
figs = []
@showprogress for gname in gnames
  f = ncp_fig(gname,(0.8,6e4))
  plot!(f,title=gname,ylabel="Conductance")
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
figs[11]
figs[12]
figs[13]
figs[14]
figs[15]
figs[16]
figs[17]
figs[18]
figs[19]



cut_ncp_fig(gnames[8],(0.8,3.5e4),(0.8,1e6))

cutfigs = []
@showprogress for gname in gnames
  f = cut_ncp_fig(gname,(0.8,3.5e4),(0.8,1e6))
  plot!(f,title=gname)
  push!(cutfigs,f)
end

plot!(figs[2],ylabel="",
  title="L10")

plot!(figs[3],ylabel="",
  title="Flickr")

plot!(figs[4],ylabel="",
  title="Commutes")

plot!(figs[5],ylabel="",
  title="Slashdot")



#study figures 
f = makefig("study-20-draft-150.smat"; title="", nlines=2,color=cgrad(cmap("CBL2")[75:230]),
          xbounds=(0.8,400))



# f = makefig("rewired-10000.0-study-20-draft-150.smat"; title="", nlines=2,color=cgrad(cmap("CBL2")[75:230]),
#           xbounds=(0.8,300))
# plot!(f,size=(400,250),dpi=1000,top_margin=-20*Measures.mm)
# Plots.savefig(f,"code/paper-figs/example-figs/epi-ncp-rewired-10000.0-study-20-draft-150.png")

# f = makefig("er-10000.0-study-20-draft-150.smat"; title="", nlines=2,color=cgrad(cmap("CBL2")[75:230]),
#           xbounds=(0.8,300))
# plot(f,size=(450,300))


#for lfr figure 

gname = "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17.smat"
f = makefig(gname;
          title="", nlines=2,color=cgrad(cmap("CBL2")[75:230]),
          xbounds=(0.8,200000),tickfontsize=18)
plot!(f,size=(900,600),dpi=1000)
Plots.savefig(f,"code/paper-figs/example-figs/epi-ncp-$(gname[1:end-5]).png")


gname = "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5.smat"
f = makefig(gname;
          title="", nlines=2,color=cgrad(cmap("CBL2")[75:230]),
          xbounds=(0.8,100000),tickfontsize=18)
plot!(f,size=(900,600),dpi=1000)
Plots.savefig(f,"code/paper-figs/example-figs/epi-ncp-$(gname[1:end-5]).png")


#for spatial figure
gname = "study-11-2022-1.smat"
f = makefig(gname;
          title="", nlines=3,color=cgrad(cmap("CBL2")[75:230]),
          xbounds=(0.8,40000),tickfontsize=18)
plot!(f,size=(600,400),dpi=1000)
Plots.savefig(f,"code/paper-figs/example-figs/epi-ncp-$(gname[1:end-5]).png")
          

gname = "study-11-2022-50.smat"
f = makefig(gname;
          title="", nlines=3,color=cgrad(cmap("CBL2")[75:230]),
          xbounds=(0.8,40000),tickfontsize=18)
plot!(f,size=(600,400),dpi=1000)
Plots.savefig(f,"code/paper-figs/example-figs/epi-ncp-$(gname[1:end-5]).png")

gname = "study-11-2023-0-noweak.smat"
f = makefig(gname;
          title="", nlines=3,color=cgrad(cmap("CBL2")[75:230]),
          xbounds=(0.8,40000),tickfontsize=18)
plot!(f,size=(600,400),dpi=1000)
Plots.savefig(f,"code/paper-figs/example-figs/epi-ncp-$(gname[1:end-5]).png")


gnames = ["study-11-2023-0-noweak.smat",
  "study-11-2022-1.smat",
  "study-11-2022-10.smat",
  "study-11-2022-20.smat",
  "study-11-2022-30.smat",
  "study-11-2022-40.smat",
  "study-11-2022-50.smat"
  ]
for gname in gnames
  f = makefig(gname;
          title="", nlines=3,color=cgrad(cmap("CBL2")[75:230]),
          xbounds=(0.8,40000),tickfontsize=18)
  plot!(f,size=(600,400),dpi=1000)
  Plots.savefig(f,"code/paper-figs/example-figs/epi-ncp-$(gname[1:end-5]).png")
end


## mexico city for missed sets figure
gname = "modmexico-city.smat"
f = makefig(gname;
          title="", nlines=4,color=cgrad(cmap("CBL2")[75:230]),
          xbounds=(0.8,2e5),ybounds = (1e-5,1.15)
          ,tickfontsize=15)
plot!(f,size=(600,400),dpi=1000,colorbar=true,
        bottom_margin=-5Measures.mm)
Plots.savefig(f,"code/paper-figs/example-figs/epi-ncp-$(gname[1:end-5]).png")


gname = "modmexico-city.smat"
f = makefig(gname;
          title="", nlines=4,color=cgrad(cmap("CBL2")[75:230]),
          xbounds=(0.8,2e5),ybounds = (1e-5,1.15)
          ,tickfontsize=15)
plot!(f,size=(600,400),dpi=1000,colorbar=true,
        bottom_margin=-5Measures.mm)




#################### FINAL VERSION ###############################
hs = ["commutes-all","mexico", "filtered", #row1
    "cn-moduillinois", "cn-Penn", "cn-modWiscon", #row2...
    "dblp","enron","anon", #row 3
    "cit-HepPh", "slashdot", "flickr", 
    "geometric","study-25-150","cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5.smat"]
    #  "rewired-10000.0-modmexico","rewired-10000.0-cn-moduillinois","er-10000.0-dblp-cc"]
gnames = map(h->h=="" ? "" : getgnames(h,"pipeline/graphs/")[1], hs) 
titles = ["US Commutes", "Mexico City\nTrace", "Filtered\nUS Flows",
    "Sparsified\nCollege-Illinois", "Sparsified\nCollege-Penn", "Sparsified\nCollege-Wisc.",    
    "Collaboration", "Email", "Facebook\nInteractions",
    "Citation", "Slashdot","Flickr",    
    "Local\nGeometric", "Geometric\nCommunities", "Random Walk\nCommunities"]
nlines=[3,3,3, 2,2,2, 2,2,2, 2,2,2, 3,3,3]

c = cgrad(cmap("CBL2")[75:230])
fs = map(x -> makefig(x[1];title=x[2],nlines=x[3],color=c), zip(gnames, ["" for i in titles],nlines))

# heights = map(nl -> 15+(nl+1)*35, [2,2,2,2,1.5]);#,0.5]);
# totalh = sum(heights)
# heights[end]=round(Int,heights[end]*1.31)
# heights = heights /sum(heights)

heights = [3.,2,2,2,3]
heights=heights./sum(heights)

annotation_locations = [
  (200,0.005), (1200,0.004), (4000,0.005) ,#row1
  (300,0.06),(300,0.06),(300,0.06), #row2
  (2000,0.03),(1000,0.03),(1000,0.08), #row3
  (200,0.05),(200,0.05),(1000,0.05),
  (100,0.03),(500,0.008),(2000,0.008)
]

ps = deepcopy(fs)
for (i,f) in enumerate(ps)
  xann,yann = annotation_locations[i]
  annotate!(f,(xann,yann,text(titles[i], 14, :center, RGB(0.,0.,0.), :top)))
end

#plotting
newf = Plots.plot(ps...,
        layout=grid(5,3,heights=heights,widths=[0.33,0.33,0.3]/(0.33+0.33+0.3)),
        size=(1000,1200),
        top_margin=-1Measures.mm,
        bottom_margin=-7Measures.mm,
        left_margin=-7Measures.mm,
        right_margin=-2Measures.mm)

plot!(newf,dpi=1000)
savefig(newf, "code/paper-figs/fig2/epidemic-ncps-final.png")
savefig(newf, "code/paper-figs/fig2/epidemic-ncps-final.pdf")

#saving individual plots 
for (ind,f) in enumerate(ps)
  gname = gnames[ind]
  ff = deepcopy(f)
  plot!(ff,size=(300,1200/5),dpi=1000)
  Plots.savefig(ff,"code/paper-figs/fig2/individual-plots/epidemic-ncp-$(gname[1:end-5]).png")
end


#other layout 
annotation_locations = [
  (1200,0.004), (1200,0.004), (4000,0.005) ,#row1
  (300,0.06),(300,0.06),(300,0.06), #row2
  (3000,0.03),(1000,0.03),(1000,0.08), #row3
  (200,0.05),(200,0.05),(1000,0.05),
  (500,0.005),(5000,0.005),(2000,0.008)
]

ps = deepcopy(fs)
for (i,f) in enumerate(ps)
  xann,yann = annotation_locations[i]
  annotate!(f,(xann,yann,text(titles[i], 14, :center, RGB(0.,0.,0.), :top)))
end

inds = [1;2;3;13;14;4;7;9;10;12]
newf = Plots.plot(ps[inds]...,
        layout=grid(2,5,heights=[0.5,0.5],widths=[1/5 for i=1:5]),
        size=(1500,500),
        bottom_margin=-5Measures.mm,
        left_margin=-9Measures.mm,
        right_margin=-2Measures.mm)

plot!(newf,dpi=300)
savefig(newf, "code/paper-figs/fig2/epidemic-ncps-v2.png")



###############rewired version of the above figure ########################
hs = ["commutes-all","mexico", "filtered", #row1
    "cn-moduillinois", "cn-Penn", "cn-modWiscon", #row2...
    "dblp","enron","anon", #row 3
    "cit-HepPh", "slashdot", "flickr", 
    "geometric","study-11-45","cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.9-0.0-100-5.smat"]
    #  "rewired-10000.0-modmexico","rewired-10000.0-cn-moduillinois","er-10000.0-dblp-cc"]
gnames = map(h->h=="" ? "" : getgnames(h,"pipeline/graphs/")[1], hs) 
gnames = map(x->"rewired-10000.0-$x",gnames) 
titles = ["US Commutes", "Mexico City\nTrace", "Filtered\nUS Flows",
    "Sparsified\nCollege-Illinois", "Sparsified\nCollege-Penn", "Sparsified\nCollege-Wisc.",    
    "Collaboration", "Email", "Facebook\nInteractions",
    "Citation", "Slashdot","Flickr",    
    "Local\nGeometric", "Geometric\nCommunities", "Random Walk\nCommunities"]
nlines=[1,1,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1]

c = cgrad(cmap("CBL2")[75:230])
fs = map(x -> makefig(x[1];title=x[2],nlines=x[3],color=c), zip(gnames, ["" for i in titles],nlines))

heights = [1,1,1,1,1]
heights = heights./sum(heights)

annotation_locations = [
  (200,0.4), (200,0.4), (200,0.4) ,#row1
  (500,0.2),(500,0.2),(500,0.2), #row2
  (200,0.3),(200,0.3),(200,0.3), #row3
  (200,0.3),(200,0.3),(400,0.3),
  (200,0.3),(200,0.3),(500,0.3)
]

ps = deepcopy(fs)
for (i,f) in enumerate(ps)
  xann,yann = annotation_locations[i]
  annotate!(f,(xann,yann,text(titles[i], 14, :center, RGB(0.,0.,0.), :top)))
end

#plotting
newf = Plots.plot(ps...,
        layout=grid(5,3,heights=heights,widths=[0.33,0.33,0.28]/(0.33+0.33+0.28)),
        size=(1000,1200),
        top_margin=-1Measures.mm,
        bottom_margin=-7Measures.mm,
        left_margin=-7Measures.mm,
        right_margin=-2Measures.mm)


plot!(newf,dpi=200)
savefig(newf, "code/paper-figs/fig2/rewired-epidemic-ncps-final.png")

############ ER rewired version ################################
hs = ["commutes-all","mexico", "filtered", #row1
    "cn-moduillinois", "cn-Penn", "cn-modWiscon", #row2...
    "dblp","enron","anon", #row 3
    "cit-HepPh", "slashdot", "flickr", 
    "geometric","study-11-45","cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.9-0.0-100-5.smat"]
    #  "rewired-10000.0-modmexico","rewired-10000.0-cn-moduillinois","er-10000.0-dblp-cc"]
gnames = map(h->h=="" ? "" : getgnames(h,"pipeline/graphs/")[1], hs) 
gnames = map(x->"er-10000.0-$x",gnames) 
titles = ["US Commutes", "Mexico City\nTrace", "Filtered\nUS Flows",
    "Sparsified\nCollege-Illinois", "Sparsified\nCollege-Penn", "Sparsified\nCollege-Wisc.",    
    "Collaboration", "Email", "Facebook\nInteractions",
    "Citation", "Slashdot","Flickr",    
    "Local\nGeometric", "Geometric\nCommunities", "Random Walk\nCommunities"]
nlines=[1,1,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1]

c = cgrad(cmap("CBL2")[75:230])
fs = map(x -> makefig(x[1];title=x[2],nlines=x[3],color=c), zip(gnames, ["" for i in titles],nlines))

# heights = map(nl -> 15+(nl+1)*35, [2,2,2,2,1.5]);#,0.5]);
# totalh = sum(heights)
# heights[end]=round(Int,heights[end]*1.31)
# heights = heights /sum(heights)

heights = [1,1,1,1,1]
heights = heights./sum(heights)

annotation_locations = [
  (200,0.4), (200,0.4), (200,0.4) ,#row1
  (200,0.2),(200,0.2),(500,0.2), #row2
  (200,0.3),(200,0.3),(200,0.3), #row3
  (200,0.3),(200,0.3),(400,0.3),
  (200,0.3),(200,0.3),(500,0.3)
]

ps = deepcopy(fs)
for (i,f) in enumerate(ps)
  xann,yann = annotation_locations[i]
  annotate!(f,(xann,yann,text(titles[i], 14, :center, RGB(0.,0.,0.), :top)))
end

#plotting
newf = Plots.plot(ps...,
        layout=grid(5,3,heights=heights,widths=[0.33,0.33,0.28]/(0.33+0.33+0.28)),
        size=(1000,1200),
        top_margin=-1Measures.mm,
        bottom_margin=-7Measures.mm,
        left_margin=-7Measures.mm,
        right_margin=-2Measures.mm)

plot!(newf,dpi=200)
savefig(newf, "code/paper-figs/fig2/er-epidemic-ncps-final.png")






#making fugres for LFR  and STUDY graphs (synthetic models fig)

gname = getgnames("study-11-2022-1","input/graphs/")[1]
makefig(gname; title="", nlines=3,color=cgrad(cmap("CBL2")))

gname = getgnames("study-11-2022-50","input/graphs/")[1]
makefig(gname; title="", nlines=3,color=cgrad(cmap("D06")))

gname = getgnames("lfr","input/graphs/")[1]
makefig(gname; title="", nlines=3,color=cgrad(cmap("D06")))

gname = "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5.smat"
makefig(gname; title="", nlines=3,color=cgrad(cmap("D06")))



##
#makegridlines!(f, [1,2,3,4,5], [-1])

## Test out bug

f = plot(margin=-20Measures.mm, topmargin = -1Measures.mm, size=(250,120), background_color_inside=:lightgrey, 
  ylims=(0.001,1), framestyle=:none, xlims=(1,1e5))
title!(f,"Test", 
  titlelocation=:left, titlevalign=:top,
  titlefontsize=12, titlefont="Avantgarde")
##
plot(
  plot(size=(250,15),title="Test",framestyle=:none,titlefontsize=12,background=:black),
  makefig("dblp-cc.smat"),
  layout=grid(2,1, heights=(0.5,0.5)),
  size=(250,165)
)

##
plot(size=(250,15),title="Test",framestyle=:none,titlefontsize=12,background=:black)

##
1+1

##
Plots.savefig(newf,)

#################################################
gname = "study-20-draft-150.smat"
x,y = load_ncp(gname)
#place ncp data on plot 
f = plot(margin=-20Measures.mm, size=(450,300), 
              ylims=(8e-3,1.5), framestyle=:none, xlims=(1,500), xscale=:log10, yscale=:log10)
  
myncpplot_conditional!(f,x,y;
    xlims=(1.3,500),ylims=(8e-3,1.5),nbins=(100,80),#
    plotmin=false, markersize=0, linecolor=:orange, linewidth=2,
    plotmedian=false,color=cgrad(cmap("D06")))

plot!(f;ylims=ylims)

f = makefig(""; 
          xbounds = (1,500),
          ybounds = (1e-3,1.2),
          title="Study-20-draft-150",
          nlines=2,
          nbins=(100,80),
          color=cgrad(cmap("D06")))

plot(f,
      size=(500,400),framestyle=:none)


hexhist = HexBinPlots.fit(HexBinPlots.HexHistogram,log10.(x),log10.(y),100,
  xlims = log10.(extrema(x)),ylims=log10.(extrema(y)))

h,vh = HexBinPlots.make_shapes(hexhist)
vmax = maximum(vh)
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
# @show vh #no need to display this
xsh
vsh
#normalize the vertical slices. 
xsh[[8*k+1:7 for k=0:length(h)-1]]
[8*k+1:7 for k=0:length(h)-1]
#Hexagons whose centroid (x-component) are the same are normalized over

f

f = Plots.plot(10.0.^xsh,10.0.^ysh,fill_z=log10.(new_vals.+1),linecolor=nothing,
      seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=false,
      xlims=xlims,ylims=ylims)
# return xsh,ysh,h,vh,f
f

xsh
vsh
new_vals









include("../graph-io.jl")

gnames = getgnames("enron","pipeline/graphs/")

A = loadGraph(gname[1],"pipeline/graphs/")
B = loadGraph(gname[end],"pipeline/graphs/")


nnz(A)
nnz(B)/nnz(A)

length(collect(triangles(A)))
length(collect(triangles(B)))






##############################################################
############### Missed Sets Companion Figures ###############
##############################################################
epi_color = cgrad(cmap("CBL2")[75:230])
f = makefig(gnames[1],nlines=nlines_vec[1],color=epi_color,
          tickfontsize=15,
          xbounds=(0.8,2e5))

plot!(f,size=(600,400),ylims=(1e-3,1.1))
fvec[1]
fvec[2]

xlims(fvec[1])


