## Script for watching evolution of embedded graph over time 

using SparseArrays
using DelimitedFiles

function graphlines(A::SparseMatrixCSC,xy::AbstractArray{T,2}) where T
  px,py = zeros(T,0),zeros(T,0)
  P = [px,py]
  rows = rowvals(A)
  skip = NaN.*xy[:,begin] # first row
  for j=1:size(A,2) # for each column
    for nzi in nzrange(A, j)
      i = rows[nzi]
      if i > j
        push!.(P, @view xy[:,i])
        push!.(P, @view xy[:,j])
        push!.(P, skip)
      end
    end
  end
  return px, py
end



gname = "study-11-2023-1-longrange-1.smat"
gname = "study-11-2023-0-noweak.smat"
A = loadGraph(gname,"input/graphs/")

#all graphs have same xy coords 
Axy = readdlm("input/graphs/study-11-2023-1-longrange-1.xy",'\t')

#epidemics 
E = EventEpidemic.SEIRData(A,beta=1.2e-2)

seed_node = rand(1:lastindex(A,1))
l,E = EventEpidemic.epidemic(E,seed_node)
sum(E.snodes)

#plottings 
px,py = graphlines(A,Axy')

plot([0;0;1;1;0],[0;1;1;0;0],leg=false,framestyle=:none,
    margins=-10Measures.mm)

f = plot(px,py, color=RGBA(0,0.66,0.8,0.5),linewidth=0.1,leg=false,
    framestyle=:none,
    margins=-10Measures.mm)

#minimize ink by only showing infected 
inds = (!).(E.snodes)
scatter!(f, Axy[:,1][inds],Axy[:,2][inds],markersize=0.4, color=RGBA(0.9,0.0,0.0,1),framestyle=:none,
    markerstrokewidth=0)


gname
tmp = "study-11-2023-1-longrange-10.smat"
ind = findfirst(eig_data[:,1].==tmp)
lam = eig_data[ind,2]

