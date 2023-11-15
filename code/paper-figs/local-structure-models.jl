#new synthetic model 

# spatial model w long range edges

#different strategies
    # random - CL model - no need to check
    # length ~ d1 -> length ~ d2


### try new synthetic model 
function spatial_graph_edges1(n::Integer,d::Integer;degreedist=LogNormal(log(4),1))
    xy = rand(d,n).- 0.5
    for col = 1:n
      xy[:,col]./=norm(xy[:,col])
    end
    T = BallTree(xy)
    # form the edges for sparse
    ei = Int[]
    ej = Int[]
    for i=1:n
      deg = min(ceil(Int,rand(degreedist)),n-1)
      idxs, dists = knn(T, xy[:,i], deg+1)
      for j in idxs
        if i != j
          push!(ei,i)
          push!(ej,j)
        end
      end
    end
    return xy, ei, ej
  end
  """
  n nodes in dimension d (random [0,1] box)
  degreedist = degree distribution function
  p = number of expected edges per node to add at ranodm (default = 0 )
  """
  function spatial_network1(n::Integer, d::Integer; degreedist=LogNormal(log(3),1),
      p::Real=0.0)
    xy, ei, ej = spatial_graph_edges1(n, d;degreedist=degreedist)
    A = sparse(ei,ej,1,n,n)
    C = A
    if p > 0
      C .+= triu(sprand(Bool, n,n,p/n),1)
      fill(C.nzval,1)
    end
    return max.(C,C'), xy
  end
  
  #keep track of edge lengths 
  
  function initial_edge_lengths(A::SparseMatrixCSC,xy)
    W = SparseMatrixCSC{Float64}(deepcopy(A))
    rowval = rowvals(A)
    for i = 1:lastindex(A,1)
      for j in rowval[nzrange(A,i)]
        d = (xy[1,i]-xy[1,j])^2 + (xy[2,i]-xy[2,j])^2
        W[i,j] = d
        W[j,i] = d
      end
    end
    return W
  end
  
  
  #rewire edges to roughly preserve length 
  # (i,j,d_ij)
  function _undirected_edges(A::SparseMatrixCSC;use_diag::Bool = false)
    rowval = rowvals(A)
    n = size(A,2)
    ndiag = nnz(diag(A))
    nedges = (nnz(A) - ndiag)/2
    if use_diag
      nedges += ndiag
    end
    edges = Vector{Tuple{Int,Int,Union{Float64,Int}}}(undef, Int(nedges))
    curedge = 0
    for j=1:n
      for nzi = nzrange(A,j)
        i = rowval[nzi]
        if j > i || (use_diag && i == j)
          curedge += 1
          edges[curedge] = (i,j,A[i,j])
        end
      end
    end
    @assert curedge == nedges
    return edges
  end
  
  #degree preserving rewiring 
  function rewire_edges_spatial(edges::Vector{Tuple{T,T}},
      k::Integer=2*ceil(Int,length(edges)*log(length(edges)))) where T <: Integer
    for rewire_step = 1:k
      eij = rand(1:length(edges))
      ers = rand(1:length(edges))
      i,j,dist_ij = edges[eij]
      r,s,dist_rs = edges[ers]
      edges[eij] = (i,s)
      edges[ers] = (r,j)
    end
    return edges
  end
  


