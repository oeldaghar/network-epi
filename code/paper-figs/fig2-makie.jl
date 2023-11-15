
##
using CairoMakie

function graphSize(gname)
    A = loadGraph(joinpath("pipeline/graphs/",gname))
    return size(A,1)
end
 
##
function load_ncp(gname,datadir="pipeline/data/")
    n = graphSize(gname)
    
    base_graph = canonical_graph_name(gname)
    
    ncp,headerinfo = readdlm(joinpath(datadir,base_graph[1:end-5],"ncpdata","ncpinfo-epidemic-seir-60000-$(gname[1:end-5]).txt"),',',header=true)
    ncp = DataFrame(ncp,vec(headerinfo))
    
    x = map(x->min(ncp.size[x],n-ncp.size[x]),1:size(ncp,1))
    y = ncp.cond
  
    return x,y
  end 

function maketestfigaxis(lg, gname; title="", nlines = 3, ybounds = (1e-3,1), xbounds = (0.8,200000), nbins=(50,30))
  
end


# make it in 2x mode... 
fig = Figure(resolution=(650*2,652*2))
g = fig[1,1] = GridLayout(7,3)
rowgap!(g, 3*2)
colgap!(g, 5*2)
for i=1:7
  for j=1:3
    #ga = g[i,j] = GridLayout(1,1)
    f = maketestfigaxis(g[i,j], "dblp"; title="Citation and Really Long")
  end
end
fig

##

function makegridlines!(ax,xs,ys,ylast;color=:lightgrey,linewidth=1)
  #jerry rig this to work for now. 
  #cant seem to find in the docs...

  snscript = ["⁻¹", "⁻²", "⁻³"]
  spscript = ["¹", "²", "³", "⁴", "⁵"]
  #vertical lines
  for x in xs
      vlines!(ax,[10.0^x];color,linewidth)
      text!(ax, 10^x, ylast; text="10$(spscript[x])", fontsize=20, color=RGB(0.75,0.75,0.75), align=(:right,:bottom))
      #annotate!(f, (10.0^x, ylims(f)[1], text("10^{$x}", 9, :right, RGB(0.75,0.75,0.75), :bottom)))
  end

  #horizontal lines
  for y in ys
      hlines!(ax,[10.0^(y)];color,linewidth)
      text!(ax, 1, 10.0^y; text="10$(snscript[-y])", fontsize=22, color=RGB(0.75,0.75,0.75), align=(:left,:top))
      #annotate!(f, (xlims(f)[1], 10.0^y,  text("10^{$y}", 9, :left, RGB(0.75,0.75,0.75), :top)))
  end
end

function myhexbinmakie!(ax,x,y;nbins=100,xlims=extrema(x),ylims=extrema(y))
  hexhist = fit(HexBinPlots.HexHistogram,log10.(x),log10.(y),nbins,
  xlims = log10.(xlims),ylims=log10.(ylims))

  h,vh = HexBinPlots.make_shapes(hexhist)
  vmax = maximum(vh)
  polys = Vector{Point2f}[] 
  crange = extrema(log10.(vh.+1))
  for k in eachindex(h)
    #push!(polys, map(z->Point2f(10.0^z[1],10.0^z[2]), zip(h[k].x, h[k].y)))
    poly!(ax, map(z->Point2f(10.0^z[1],10.0^z[2]), zip(h[k].x, h[k].y));
       color=[log10.(vh[k]+1)], colorrange=crange, colormap=:inferno)
  end
  
  #poly!.(ax,polys,color=log10.(vh.+1))
end

function myhexbinmakie!(ax,x,y;nbins=100,xlims=extrema(x),ylims=extrema(y))
    hexhist = fit(HexBinPlots.HexHistogram,log10.(x),log10.(y),nbins,
    xlims = log10.(xlims),ylims=log10.(ylims))
  
    h,vh = HexBinPlots.make_shapes(hexhist)
    vmax = maximum(vh)
    polys = Vector{Point2f}[] 
    crange = extrema(log10.(vh.+1))
    for k in eachindex(h)
      #push!(polys, map(z->Point2f(10.0^z[1],10.0^z[2]), zip(h[k].x, h[k].y)))
      poly!(ax, map(z->Point2f(10.0^z[1],10.0^z[2]), zip(h[k].x, h[k].y));
         color=[log10.(vh[k]+1)], colorrange=crange, colormap=:inferno)
    end
    
    #poly!.(ax,polys,color=log10.(vh.+1))
  end
  





function makefigaxis(lg, gname; title="", 
    nlines = 3, ybounds = (1e-3,1), xbounds = (0.8,200000), nbins=(50,30))
  ydim = 52*(nlines+1) + (isempty(title) ? 0 : 16*2) + (nlines == 2 ? 1 : 0)
  ylims = (10.0^(-nlines-1),1.3)
  a = Axis(lg[1,1],
    spinewidth=0,
    title=title,
    titlesize=28,
    titlealign=:left,
    titlefont=:regular,
    titlegap=0,
    xautolimitmargin=(0,0),
    yautolimitmargin=(0,0),
    xscale=log10,yscale=log10)
  CairoMakie.xlims!(a,(1,2e5))
  CairoMakie.ylims!(a,ylims)
  hidedecorations!(a)
  hidespines!(a)
  
  if gname != "" 
    makegridlines!(a, [1,2,3,4,5], -(1:nlines), ylims[1])

    x,y = load_ncp(gname)
    myhexbinmakie!(a, x, y; nbins, xlims=xbounds, ylims=ybounds )
  end
  return lg     
end

# make it in 2x mode... 

hs = ["dblp" "enron" "" #row 1
    "anon" "" "" #row2
    "cn-moduillinois" "cn-Penn" "cn-modWiscon" #row3...
    "commutes-all" "mexico" ""
    "geometric" "study-11-45" "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.9-0.0-100-5.smat"
     "rewired-10000.0-modmexico" "rewired-10000.0-cn-moduillinois" "er-10000.0-dblp-cc"]
gnames = map(h->h=="" ? "" : getgnames(h,"pipeline/graphs/")[1], hs) 
titles = ["Citation" "Email" ""
 "Facebook Interactions"  "" ""
 "Sparsified College-Illinois"  "Sparsified College-Penn" "Sparsified College-Wisconsin"
 "US Commutes" "Mexico City Trace" ""
 "Local Geometric" "Geometric Communities" "Random Walk Communities"
 "Config. Model of Mexico City Trace" "Config. Model of  Sparse. FB-Illinois" "Randomized Citation"]

nlinesperrow=[3,2,2,3,3,1] 

fig = Figure(resolution=(700*2,752*2))
g = fig[1,1] = GridLayout(7,3)
rowgap!(g, 10*2)
colgap!(g, 5*2)
for i=1:6
  for j=1:3
    #ga = g[i,j] = GridLayout(1,1)
    f = makefigaxis(g[i,j], gnames[i,j]; title=titles[i,j], nlines=nlinesperrow[i])
  end
  rowsize!(g, i, Relative(nlinesperrow[i]/sum(nlinesperrow)))
end
fig

##
save("epideic-ncps-fig-2-makie.svg", fig, pt_per_unit=0.33)