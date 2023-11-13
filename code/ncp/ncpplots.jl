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
include("hexbins.jl")
using StatsBase, Plots
function myhexbin(x,y;nbins=100)
  hexhist = fit(HexBinPlots.HexHistogram,log10.(x),log10.(y),nbins)
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
  @show vh
  Plots.plot(10.0.^xsh,10.0.^ysh,fill_z=log10.(vh.+1),linecolor=nothing,
        seriestype=:shape,xscale=:log10,yscale=:log10,label="",colorbar=false)
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

  cv = CategoricalArrays.cut(log10.(x),floor(Int,nbins/4);allowempty=true)
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
