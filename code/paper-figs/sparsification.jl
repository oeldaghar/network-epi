## for generating a figure showing impact of sparsification

mainDir = "/p/mnt/scratch/network-epi/code/"
include(joinpath(mainDir,"graph-io.jl"))
include(joinpath(mainDir,"data-io.jl"))

include(joinpath(mainDir,"ncp/ncp-acl.jl"))

using Plots
using Measures

#original graph
gname = "Penn94.smat"
A = loadGraph(gname,"input/graphs/")

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

dst = "/p/mnt/scratch/network-epi/pipeline/data/Penn94/ncpdata/"
ncp,headerinfo = readdlm(dst*"ncpinfo-$(gname[1:end-5]).txt",',',header=true)
ncp = DataFrame(ncp,vec(headerinfo))

x,y = ncp.size,ncp.cond
nnodes = size(A,1)

x = map(a->min(a,nnodes-a),x)
xbounds,ybounds = (1.0,5e4),(1e-3,1.15)

f = plot(margin=-20Measures.mm,  topmargin=isempty("") ? 0Measures.mm : -nlines*0.5Measures.mm , size=(600,400), 
                ylims=ybounds, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
makegridlines!(f, [1,2,3,4], -(1:2),tickfontsize=15)

c = cgrad(cgrad(:viridis)[1:250])

myhexbin2!(f,x,y,xlims=xbounds,ylims=ybounds,color=c,nbins=80,normalize_color=true)
myhexbin2!(f,x,y,xlims=xbounds,ylims=ybounds,color=c,nbins=80,normalize_color=true)

plot!(f,colorbar=true,
        bottom_margin=-5Measures.mm,
        right_margin=-1Measures.mm,
        dpi=1000)
savefig(f,"/p/mnt/scratch/network-epi/code/paper-figs/example-figs/acl-ncp-$(gname[1:end-5]).png")



## sparsified graph
gname = "cn-Penn94.smat"
A = loadGraph(gname,"input/graphs/")

dst = "/p/mnt/scratch/network-epi/pipeline/data/cn-Penn94/ncpdata/"
ncp,headerinfo = readdlm(dst*"ncpinfo-$(gname[1:end-5]).txt",',',header=true)
ncp = DataFrame(ncp,vec(headerinfo))

x,y = ncp.size,ncp.cond
nnodes = size(A,1)

x = map(a->min(a,nnodes-a),x)

f = plot(margin=-20Measures.mm,  topmargin=isempty("") ? 0Measures.mm : -nlines*0.5Measures.mm , size=(600,400), 
                ylims=ybounds, framestyle=:none, xlims=xbounds, xscale=:log10, yscale=:log10)
makegridlines!(f, [1,2,3,4], -(1:2),tickfontsize=15)

c = cgrad(cgrad(:viridis)[1:250])

myhexbin2!(f,x,y,xlims=xbounds,ylims=ybounds,color=c,nbins=80,normalize_color=true)
myhexbin2!(f,x,y,xlims=xbounds,ylims=ybounds,color=c,nbins=80,normalize_color=true)

plot!(f,colorbar=true,
        bottom_margin=-5Measures.mm,
        right_margin=-1Measures.mm,
        dpi=1000)
savefig(f,"/p/mnt/scratch/network-epi/code/paper-figs/example-figs/acl-ncp-$(gname[1:end-5]).png")




