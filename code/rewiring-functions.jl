using MatrixNetworks, LinearAlgebra, Random, SparseArrays#, Sparsification


"""
    _undirected_edges(A::SparseMatrixCSC;use_diag::Bool = false)

given a graph A, retrieve undirected edges from A for rewiring
"""
function _undirected_edges(A::SparseMatrixCSC;use_diag::Bool = false)
  rowval = rowvals(A)
  n = size(A,2)
  ndiag = nnz(diag(A))
  nedges = (nnz(A) - ndiag)/2
  if use_diag
    nedges += ndiag
  end
  edges = Vector{Tuple{Int,Int}}(undef, Int(nedges))
  curedge = 0
  for j=1:n
    for nzi = nzrange(A,j)
      i = rowval[nzi]
      if j > i || (use_diag && i == j)
        curedge += 1
        edges[curedge] = (i,j)
      end
    end
  end
  @assert curedge == nedges
  return edges
end

#degree preserving rewiring 
function rewire_edges(edges::Vector{Tuple{T,T}},
    k::Integer=2*ceil(Int,length(edges)*log(length(edges)))) where T <: Integer
  for rewire_step = 1:k
    eij = rand(1:length(edges))
    ers = rand(1:length(edges))
    i,j = edges[eij]
    r,s = edges[ers]
    edges[eij] = (i,s)
    edges[ers] = (r,j)
  end
  return edges
end

function rewire_graph(A::SparseMatrixCSC, k::Integer=ceil(Int,nnz(A)*log(nnz(A))))
  # @assert is_undirected(A)
  new_edges = rewire_edges(_undirected_edges(A),k)
  B = sparse(map(first, new_edges), map(x->x[2], new_edges), 1, size(A)...)
  fill!(B.nzval,1.0)
  return max.(B,B')
end

#er rewiring 
function er_rewire_edges(nnodes::Int,edges::Vector{Tuple{T,T}},
  k::Integer=2*ceil(Int,length(edges)*log(length(edges)))) where T <: Integer
  """erdos renyi rewiring, sample an edge and randomly choose its endpoints"""
  for rewire_step = 1:k
    u,v = rand(1:nnodes,2)
    while u==v
      u,v = rand(1:nnodes,2)
    end
    edges[rand(1:length(edges))] = (u,v) #tuple(rand(1:nnodes,2)...)
  end
  return edges
end

function er_rewire_graph(A::SparseMatrixCSC, k::Integer=ceil(Int,nnz(A)*log(nnz(A))))
  # @assert is_undirected(A)
  new_edges = er_rewire_edges(size(A,1),_undirected_edges(A),k)
  B = sparse(map(first, new_edges), map(x->x[2], new_edges), 1, size(A)...)
  dropzeros!(B)
  fill!(B.nzval,1.0)
  return max.(B,B')
end


#triangle rewiring 
#not quite faithful to triangles and degree after projecting to pariwise
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

