# rough triangle preserving rewiring. not quite what we want but it'll do.

mainDir = joinpath("/",split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))]...)
include(joinpath(mainDir,"code/graph-io.jl"))

using SparseMatrixDicts, MatrixNetworks, SparseArrays, ProgressMeter,Random 
function random_triplet(n::Int)
  while true
    i = rand(1:n)
    j = rand(1:n)
    k = rand(1:n)
    if i != j != k
      return i,j,k
    end
  end
  return 0,0,0
end
function add_triangle!(W::SparseMatrixDict,t)
  u,v,w = t
  W[u,v] += 1
  W[v,u] += 1
  W[u,w] += 1
  W[w,u] += 1
  W[v,w] += 1
  W[w,v] += 1
end
function remove_triangle!(W::SparseMatrixDict,t)
  u,v,w = t
  W[u,v] -= 1
  W[v,u] -= 1
  W[u,w] -= 1
  W[w,u] -= 1
  W[v,w] -= 1
  W[w,v] -= 1
end
function swap_triangles(A::SparseMatrixCSC, nsteps::Int)
  # build the list of triangles...
  # A = sparse(first.(edges), last.(edges), 1, n, n)
  # A = max.(A,A')
  n = size(A,1)
  W = SparseMatrixDict(n,n)
  tris = collect(triangles(A))
  for tri in tris
    u,v,w = tri
    add_triangle!(W, (u,v,w))
  end
  B = sparse(W)
  B = min.(B,1)
  C = dropzeros!(A .- B) # these are all the edges not in a triangle
  for (i,j,v) in zip(findnz(C)...)
    W[i,j] += v
  end

  # now swap triangles...
  for step=1:nsteps
    t1id = rand(1:length(tris))
    t2 = random_triplet(n)
    t = tris[t1id]
    remove_triangle!(W,t)
    add_triangle!(W,t2)
    tris[t1id] = t2
  end
  R = sparse(W)
  R = min.(R,1) # remove weights
  return R
end

# gs = ["flickr","livejournal"]
# total_steps = 1000000
# rseed = 7
# Random.seed!(rseed)
# #ps = collect(0.1:0.1:1.0) 
# ps = [1.0] #for flickr  + livejournal
# seeds = abs.(rand(Int,length(ps)))
# gpath = "../data/graphs/"
# dst = gpath
# for g in gs
#   println("working on graph $g")
#   gname = getgnames(g,gpath)[1]
#   A = loadGraph(gname,gpath)
#   @showprogress for (i,p) in enumerate(ps)
#     Random.seed!(seeds[i])
#     B = swap_triangles(A,Int(p*total_steps))
#     writeSMAT(B,dst*"triangle-rewired-$p-$rseed-"*gname)
#   end
# end
