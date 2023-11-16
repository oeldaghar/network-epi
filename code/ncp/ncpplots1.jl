#=
the code below is mostly attributed to David. I've only made slight modifications.
=#


#=
NCPPlots

Goal, make some nicer looking NCP plots like we had in the LocalGraphClustering repo.

An NCP plot in the LGC repo had:

- hexbin for histogram2d with log-spaced x, y values
- but real labels.

(Copyright note, I did these all from memory, and did not check the
LGC code...)

https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.cut.html
=#

# https://stackoverflow.com/questions/47391489/group-dataframe-by-binning-a-columnfloat64-in-julia
# Also Categorical Arrays cuts
##
using CategoricalArrays
# cut(log10.(ncpdata.cond), 50)
using Plots

##
#allows for execution from command line as well as an ide so files can be run modularlly
mainDir = joinpath("/",split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))]...)

include(joinpath(mainDir,"code","ncp","hexbins1.jl"))

using StatsBase, Plots
function myhexbin(x,y;nbins=100)
  hexhist = HexBinPlots.fit(HexBinPlots.HexHistogram,log10.(x),log10.(y),nbins)
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
  f = Plots.plot(10.0.^xsh,10.0.^ysh,fill_z=log10.(vsh.+1),linecolor=nothing,
        seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=false)
  f
end

function myhexbin_local(x,y;nbins=100)
  hexhist = HexBinPlots.fit(HexBinPlots.HexHistogram,log10.(x),log10.(y),nbins)
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

  #using fill_z appears to do something different locally vs on the server.
  #so plot each bin one-by-one.
  f = Plots.plot(10.0.^xsh[1:8],10.0.^ysh[1:8],fill_z=log10(vsh[1]+1),linecolor=nothing,
        seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=false,
        ylims=(1e-4,1),xlims=(3,2*1.2*maximum(x)),leg=false)
  f
end
# x=[1,2.0]
# y=[1,2.0]
# myhexbin(x,y,nbins=50)
# scatter!(x,y,legend=false)


##
#myhexbin(randn(1000).^2,randn(1000).^2)

using CategoricalArrays, Statistics
# regarding groupby
# https://github.com/JuliaLang/julia/issues/32331
function myncpplot(x,y;nbins=100,plotmin::Bool=true,plotmedian::Bool=true,xscale=:log10,yscale=:log10,
                    linestyle=:solid,markersize=4,markershape=:circle)
  f = myhexbin(x,y,nbins=nbins)
  #scatter!(x,y,alpha=0.25)

  minx = Vector{Float64}()
  miny = Vector{Float64}()
  medianx = Vector{Float64}()
  mediany = Vector{Float64}()

  cv = CategoricalArrays.cut(log10.(x),floor(Int,first(nbins)/10);allowempty=true)
  p = sortperm(cv)
  firstindex = 1
  while firstindex <= length(p)
    first = cv[p[firstindex]]
    lastindex = firstindex + 1
    while lastindex <= length(p) && cv[p[lastindex]] == first
      lastindex += 1
    end
    # get the index of the minimizing element of y
    imin = p[firstindex + argmin(@view y[p[firstindex:lastindex-1]]) - 1]
    #println(first, " ", firstindex, " ", lastindex, " ", imin)
    push!(minx, x[imin])
    push!(miny, y[imin])

    push!(medianx, median(@view x[p[firstindex:lastindex-1]]))
    push!(mediany, median(@view y[p[firstindex:lastindex-1]]))
    firstindex = lastindex # setup for next
  end
  if plotmin
    Plots.plot!(f,minx,miny,label="", color=1,xscale=xscale,yscale=yscale,
                linestyle=linestyle,markersize=markersize,markershape=markershape)
  end
  if plotmedian
    Plots.plot!(f,medianx,mediany,label="", color=1,yscale=yscale,
                linestyle=linestyle,markersize=markersize,markershape=markershape)
  end
  return f
end

function myncpplot_local(x,y;nbins=100,plotmin::Bool=true,plotmedian::Bool=true,xscale=:log10,yscale=:log10,
                    linestyle=:solid,markersize=4,markershape=:circle)
  f = myhexbin_local(x,y,nbins=nbins)
  #scatter!(x,y,alpha=0.25)

  minx = Vector{Float64}()
  miny = Vector{Float64}()
  medianx = Vector{Float64}()
  mediany = Vector{Float64}()

  cv = CategoricalArrays.cut(log10.(x),floor(Int,first(nbins)/10);allowempty=true)
  p = sortperm(cv)
  firstindex = 1
  while firstindex <= length(p)
    first = cv[p[firstindex]]
    lastindex = firstindex + 1
    while lastindex <= length(p) && cv[p[lastindex]] == first
      lastindex += 1
    end
    # get the index of the minimizing element of y
    imin = p[firstindex + argmin(@view y[p[firstindex:lastindex-1]]) - 1]
    #println(first, " ", firstindex, " ", lastindex, " ", imin)
    push!(minx, x[imin])
    push!(miny, y[imin])

    push!(medianx, median(@view x[p[firstindex:lastindex-1]]))
    push!(mediany, median(@view y[p[firstindex:lastindex-1]]))
    firstindex = lastindex # setup for next
  end
  if plotmin
    Plots.plot!(f,minx,miny,label="", color=1,xscale=xscale,yscale=yscale,
                linestyle=linestyle,markersize=markersize,markershape=markershape)
  end
  if plotmedian
    Plots.plot!(f,medianx,mediany,label="", color=1,yscale=yscale,
                linestyle=linestyle,markersize=markersize,markershape=markershape)
  end
  f
end


##test 
function myhexbin1(x,y;nbins=100,xlims=extrema(x),ylims=extrema(y))
  hexhist = fit(HexBinPlots.HexHistogram,log10.(x),log10.(y),nbins,
    xlims = log10.(xlims),ylims=log10.(ylims))

  h,vh = HexBinPlots.make_shapes(hexhist)
  vmax = maximum(vh)
  xsh = Vector{Float64}()
  ysh = Vector{Float64}()
  for k in eachindex(h)
    append!(xsh,h[k].x)
    push!(xsh,h[k].x[1])
    push!(xsh,NaN)
    append!(ysh,h[k].y)
    push!(ysh,h[k].y[1])
    push!(ysh,NaN)
  end
  # @show vh #no need to display this
  f = Plots.plot(10.0.^xsh,10.0.^ysh,fill_z=log10.(vh.+1),linecolor=nothing,
        seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=false,
        xlims=xlims,ylims=ylims)
  f
end
function myncpplot1(x,y;nbins=100,plotmin::Bool=true,plotmedian::Bool=true,xscale=:log10,yscale=:log10,
                    linestyle=:solid,markersize=4,markershape=:circle,
                    xlims=extrema(x),ylims=extrema(y))
  f = myhexbin1(x,y,nbins=nbins,xlims=xlims,ylims=ylims)
  #scatter!(x,y,alpha=0.25)

  minx = Vector{Float64}()
  miny = Vector{Float64}()
  medianx = Vector{Float64}()
  mediany = Vector{Float64}()

  cv = CategoricalArrays.cut(log10.(x),floor(Int,first(nbins)/10);allowempty=true)
  p = sortperm(cv)
  firstindex = 1
  while firstindex <= length(p)
    first = cv[p[firstindex]]
    lastindex = firstindex + 1
    while lastindex <= length(p) && cv[p[lastindex]] == first
      lastindex += 1
    end
    # get the index of the minimizing element of y
    imin = p[firstindex + argmin(@view y[p[firstindex:lastindex-1]]) - 1]
    #println(first, " ", firstindex, " ", lastindex, " ", imin)
    push!(minx, x[imin])
    push!(miny, y[imin])

    push!(medianx, median(@view x[p[firstindex:lastindex-1]]))
    push!(mediany, median(@view y[p[firstindex:lastindex-1]]))
    firstindex = lastindex # setup for next
  end
  if plotmin
    Plots.plot!(f,minx,miny,label="", color=1,xscale=xscale,yscale=yscale,
                linestyle=linestyle,markersize=markersize,markershape=markershape)
  end
  if plotmedian
    Plots.plot!(f,medianx,mediany,label="", color=1,yscale=yscale,
                linestyle=linestyle,markersize=markersize,markershape=markershape)
  end
  return f
end


###### CODE for V1.5.1
# function myhexbin_conditional(x,y;nbins=100,xlims=extrema(x),ylims=extrema(y))
#   hexhist = HexBinPlots.fit(HexBinPlots.HexHistogram,log10.(x),log10.(y),nbins,
#     xlims = log10.(xlims),ylims=log10.(ylims))

#   h,vh = HexBinPlots.make_shapes(hexhist)
#   vmax = maximum(vh)
#   xsh = Vector{Float64}()
#   ysh = Vector{Float64}()
#   for k in eachindex(h)
#     append!(xsh,h[k].x)
#     push!(xsh,h[k].x[1])
#     push!(xsh,NaN)
#     append!(ysh,h[k].y)
#     push!(ysh,h[k].y[1])
#     push!(ysh,NaN)
#   end
#   # @show vh #no need to display this

#   #normalize the vertical slices. 
#   #Hexagons whose centroid (x-component) are the same are normalized over
#   tmp = Dict()
#   for k in eachindex(h)
#     centroid = mean(h[k].x) #x-component only
#     if centroid in keys(tmp)
#       push!(tmp[centroid],k)
#     else
#       tmp[centroid] = Vector{Int}([k])
#     end
#   end

#   new_vals = zeros(Float64,lastindex(vh))
#   for k in keys(tmp)
#     inds = tmp[k]
#     new_vals[inds] .= vh[inds]/maximum(vh[inds]) #normalize by maximum
#   end

#   f = Plots.plot(10.0.^xsh,10.0.^ysh,fill_z=log10.(new_vals.+1),linecolor=nothing,
#         seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=false,
#         xlims=xlims,ylims=ylims)
#   # return xsh,ysh,h,vh,f
#   f
# end

# function myhexbin_conditional!(f,x,y;nbins=100,xlims=extrema(x),ylims=extrema(y))
#   hexhist = HexBinPlots.fit(HexBinPlots.HexHistogram,log10.(x),log10.(y),nbins,
#     xlims = log10.(xlims),ylims=log10.(ylims))

#   h,vh = HexBinPlots.make_shapes(hexhist)
#   vmax = maximum(vh)
#   xsh = Vector{Float64}()
#   ysh = Vector{Float64}()
#   for k in eachindex(h)
#     append!(xsh,h[k].x)
#     push!(xsh,h[k].x[1])
#     push!(xsh,NaN)
#     append!(ysh,h[k].y)
#     push!(ysh,h[k].y[1])
#     push!(ysh,NaN)
#   end
#   # @show vh #no need to display this

#   #normalize the vertical slices. 
#   #Hexagons whose centroid (x-component) are the same are normalized over
#   tmp = Dict()
#   for k in eachindex(h)
#     centroid = mean(h[k].x) #x-component only
#     if centroid in keys(tmp)
#       push!(tmp[centroid],k)
#     else
#       tmp[centroid] = Vector{Int}([k])
#     end
#   end

#   new_vals = zeros(Float64,lastindex(vh))
#   for k in keys(tmp)
#     inds = tmp[k]
#     new_vals[inds] .= vh[inds]/maximum(vh[inds]) #normalize by maximum
#   end

#   f = Plots.plot!(f,10.0.^xsh,10.0.^ysh,fill_z=log10.(new_vals.+1),linecolor=nothing,
#         seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=false,
#         xlims=xlims,ylims=ylims)
#   # return xsh,ysh,h,vh,f
#   f
# end

# function myncpplot_conditional!(f,x,y;nbins=100,plotmin::Bool=true,plotmedian::Bool=true,xscale=:log10,yscale=:log10,
#                     linestyle=:solid,markersize=4,markershape=:circle,linecolor=1,linewidth=2,
#                     xlims=extrema(x),ylims=extrema(y))
#   myhexbin_conditional!(f,x,y,nbins=nbins,xlims=xlims,ylims=ylims)
#   #scatter!(x,y,alpha=0.25)

#   minx = Vector{Float64}()
#   miny = Vector{Float64}()
#   medianx = Vector{Float64}()
#   mediany = Vector{Float64}()

#   cv = CategoricalArrays.cut(log10.(x),log10.([1,10,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000]);allowempty=true)
#   p = sortperm(cv)
#   firstindex = 1
#   while firstindex <= length(p)
#     first = cv[p[firstindex]]
#     lastindex = firstindex + 1
#     while lastindex <= length(p) && cv[p[lastindex]] == first
#       lastindex += 1
#     end
#     # get the index of the minimizing element of y
#     imin = p[firstindex + argmin(@view y[p[firstindex:lastindex-1]]) - 1]
#     #println(first, " ", firstindex, " ", lastindex, " ", imin)
#     push!(minx, x[imin])
#     push!(miny, y[imin])

#     push!(medianx, minimum(@view x[p[firstindex:lastindex-1]]))
#     push!(mediany, median(@view y[p[firstindex:lastindex-1]]))
#     firstindex = lastindex # setup for next
#   end
#   if plotmin
#     Plots.plot!(f,minx,miny,label="", color=linecolor,xscale=xscale,yscale=yscale,
#                 linestyle=linestyle,markersize=markersize,markershape=markershape,linewidth=linewidth)
#   end
#   if plotmedian
#     @show medianx, mediany
#     Plots.plot!(f,medianx,mediany,label="", color=linecolor,yscale=yscale,
#                 linestyle=linestyle,markersize=markersize,markershape=markershape,linewidth=linewidth)
#   end
#   return f
# end
###### END of CODE for V1.5.1

## V1.8.0
function myhexbin_conditional(x,y;nbins=100,xlims=extrema(x),ylims=extrema(y),color=:inferno,colorbar::Bool=false)
  hexhist = HexBinPlots.fit(HexBinPlots.HexHistogram,log10.(x),log10.(y),nbins,
    xlims = log10.(xlims),ylims=log10.(ylims))

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
  
  #normalize the vertical slices. 
  #Hexagons whose centroid (x-component) are the same are normalized over
  tmp = Dict()
  for k in eachindex(h)
    centroid = mean(h[k].x) #x-component only
    if centroid in keys(tmp)
      append!(tmp[centroid],collect(8*(k-1)+1:8*(k-1)+8))
    else
      tmp[centroid] = Vector{Int}(collect(8*(k-1)+1:8*(k-1)+8))
    end
  end

  new_vals = zeros(Float64,lastindex(vsh))
  for k in keys(tmp)
    inds = tmp[k]
    new_vals[inds] .= vsh[inds]/maximum(vsh[inds]) #normalize by maximum
  end

  f = Plots.plot(10.0.^xsh,10.0.^ysh,fill_z=log10.(new_vals.+1)/log10(2),linecolor=nothing,
        seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=colorbar,
        xlims=xlims,ylims=ylims,c=color)
  # return xsh,ysh,h,vh,f
  f
end

function myhexbin_conditional!(f,x,y;nbins=100,xlims=extrema(x),ylims=extrema(y),color=:inferno,colorbar=false)
  hexhist = HexBinPlots.fit(HexBinPlots.HexHistogram,log10.(x),log10.(y),nbins,
    xlims = log10.(xlims),ylims=log10.(ylims))

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

  #normalize the vertical slices. 
  #Hexagons whose centroid (x-component) are the same are normalized over
  tmp = Dict()
  for k in eachindex(h)
    centroid = mean(h[k].x) #x-component only
    if centroid in keys(tmp)
      append!(tmp[centroid],collect(8*(k-1)+1:8*(k-1)+8))
    else
      tmp[centroid] = Vector{Int}(collect(8*(k-1)+1:8*(k-1)+8))
    end
  end

  new_vals = zeros(Float64,lastindex(vsh))
  for k in keys(tmp)
    inds = tmp[k]
    new_vals[inds] .= vsh[inds]/maximum(vsh[inds]) #normalize by maximum
  end

  f = Plots.plot!(f,10.0.^xsh,10.0.^ysh,fill_z=log10.(new_vals.+1)/log10(2),linecolor=nothing,
        seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=colorbar,
        xlims=xlims,ylims=ylims,c=color)
  # return xsh,ysh,h,vh,f
  f
end

function myncpplot_conditional!(f,x,y;nbins=100,plotmin::Bool=true,plotmedian::Bool=true,xscale=:log10,yscale=:log10,
                    linestyle=:solid,markersize=4,markershape=:circle,linecolor=1,linewidth=2,
                    xlims=extrema(x),ylims=extrema(y),color=:inferno,colorbar=false)
  myhexbin_conditional!(f,x,y,nbins=nbins,xlims=xlims,ylims=ylims,color=color,colorbar=colorbar)
  #scatter!(x,y,alpha=0.25)

  minx = Vector{Float64}()
  miny = Vector{Float64}()
  medianx = Vector{Float64}()
  mediany = Vector{Float64}()

  cv = CategoricalArrays.cut(log10.(x),log10.([1,10,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000]);allowempty=true)
  p = sortperm(cv)
  firstindex = 1
  while firstindex <= length(p)
    first = cv[p[firstindex]]
    lastindex = firstindex + 1
    while lastindex <= length(p) && cv[p[lastindex]] == first
      lastindex += 1
    end
    # get the index of the minimizing element of y
    imin = p[firstindex + argmin(@view y[p[firstindex:lastindex-1]]) - 1]
    #println(first, " ", firstindex, " ", lastindex, " ", imin)
    push!(minx, x[imin])
    push!(miny, y[imin])

    push!(medianx, minimum(@view x[p[firstindex:lastindex-1]]))
    push!(mediany, median(@view y[p[firstindex:lastindex-1]]))
    firstindex = lastindex # setup for next
  end
  if plotmin
    Plots.plot!(f,minx,miny,label="", color=linecolor,xscale=xscale,yscale=yscale,
                linestyle=linestyle,markersize=markersize,markershape=markershape,linewidth=linewidth)
  end
  if plotmedian
    @show medianx, mediany
    Plots.plot!(f,medianx,mediany,label="", color=linecolor,yscale=yscale,
                linestyle=linestyle,markersize=markersize,markershape=markershape,linewidth=linewidth)
  end
  return f
end
##



# updating functions 
function myhexbin2!(f,x,y;nbins=100,xlims=extrema(x),ylims=extrema(y),color=:inferno,
                      normalize_color::Bool=false)
  
  hexhist = fit(HexBinPlots.HexHistogram,log10.(x),log10.(y),nbins,
    xlims = log10.(xlims),ylims=log10.(ylims))

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
  
  color_weight = log10.(vsh.+1)

  if normalize_color #rescale colors from (min,max) to (0,1)
    cmin,cmax = extrema(color_weight)
    crange = cmax-cmin
    color_weight = (color_weight.-cmin)./crange
  end

  Plots.plot!(f,10.0.^xsh,10.0.^ysh,fill_z=color_weight,linecolor=nothing,
        seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=false,
        xlims=xlims,ylims=ylims,color=color)
  f
end

function myncpplot2!(f,x,y;nbins=100,plotmin::Bool=true,plotmedian::Bool=true,xscale=:log10,yscale=:log10,
                    linestyle=:solid,markersize=4,markershape=:circle,linecolor=1,linewidth=2,
                    xlims=extrema(x),ylims=extrema(y))
  myhexbin2!(f,x,y,nbins=nbins,xlims=xlims,ylims=ylims)
  #scatter!(x,y,alpha=0.25)

  minx = Vector{Float64}()
  miny = Vector{Float64}()
  medianx = Vector{Float64}()
  mediany = Vector{Float64}()

  cv = CategoricalArrays.cut(log10.(x),log10.([1,10,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000]);allowempty=true)
  p = sortperm(cv)
  firstindex = 1
  while firstindex <= length(p)
    first = cv[p[firstindex]]
    lastindex = firstindex + 1
    while lastindex <= length(p) && cv[p[lastindex]] == first
      lastindex += 1
    end
    # get the index of the minimizing element of y
    imin = p[firstindex + argmin(@view y[p[firstindex:lastindex-1]]) - 1]
    #println(first, " ", firstindex, " ", lastindex, " ", imin)
    push!(minx, x[imin])
    push!(miny, y[imin])

    push!(medianx, minimum(@view x[p[firstindex:lastindex-1]]))
    push!(mediany, median(@view y[p[firstindex:lastindex-1]]))
    firstindex = lastindex # setup for next
  end
  if plotmin
    Plots.plot!(f,minx,miny,label="", color=linecolor,xscale=xscale,yscale=yscale,
                linestyle=linestyle,markersize=markersize,markershape=markershape,linewidth=linewidth)
  end
  if plotmedian
    @show medianx, mediany
    Plots.plot!(f,medianx,mediany,label="", color=linecolor,yscale=yscale,
                linestyle=linestyle,markersize=markersize,markershape=markershape,linewidth=linewidth)
  end
  return f
end

# function for finding difference in two different ncpplots1