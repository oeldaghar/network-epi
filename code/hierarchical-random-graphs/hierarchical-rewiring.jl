#=
rewire edges based on group membership
main idea is to use edge swaps that preserve node degrees but 
cache edge group membership and sample using those bins.

This is basically a partition of edges and rewiring within each bin. 
=#
using Random
using SparseArrays
using ProgressMeter


# include("nested-sbm.jl")

"""
    _get_binned_edges(A::SparseMatrixCSC,node_membership::Vector{Int})

use node_membership to bin the edges of A. this returns a dictionary of edge keys to edges
"""
function _get_binned_edges(A::SparseMatrixCSC,node_membership::Vector{T}) where T
    @assert(lastindex(node_membership)==lastindex(A,1)==lastindex(A,2),"number of nodes in A must match with node_membership")
    #convert to dict for edge format
    rows = rowvals(A)
    E_dict = Dict{Tuple}{Set}()
    @showprogress for col = 1:lastindex(A,1)
        col_key = node_membership[col]
        for row in rows[nzrange(A,col)]
            edge = (row,col)
            row_key = node_membership[row]
            edge_key = (row_key,col_key)
            if haskey(E_dict,edge_key)
                push!(E_dict[edge_key],edge)
            else
                E_dict[edge_key] = Set([edge])
            end
        end
    end
    return E_dict
end

"""
    _get_binned_edges_undirected(A::SparseMatrixCSC,node_membership::Vector{Int})

use node_membership to bin the edges of A. this returns a dictionary of edge keys to edges
"""
function _get_binned_edges_undirected(A::SparseMatrixCSC,node_membership::Vector{T}) where T
    @assert(lastindex(node_membership)==lastindex(A,1)==lastindex(A,2),"number of nodes in A must match with node_membership")
    #convert to dict for edge format
    rows = rowvals(A)
    E_dict = Dict{Tuple}{Set}()
    for col = 1:lastindex(A,1)
        col_key = node_membership[col]
        for row in rows[nzrange(A,col)]
            if row<col
                edge = (row,col)
                row_key = node_membership[row]
                edge_key = (row_key,col_key)
                if haskey(E_dict,edge_key)
                    push!(E_dict[edge_key],edge)
                else
                    E_dict[edge_key] = Set([edge])
                end
            end
        end
    end
    return E_dict
end


"""
    edges_to_graph(E_dict::Dict{Tuple}{Set})

given a dictionary of binned edges, convert this to a sparse, symmetric matrix
"""
function edges_to_graph(E_dict::Dict{Tuple}{Set})
    ei,ej = [],[]
    for key in keys(E_dict)
        bin_edges = E_dict[key]
        for edge in bin_edges
            push!(ei,edge[1])
            push!(ej,edge[2])
        end
    end
    nnodes = max(maximum(ei),maximum(ej))
    X = sparse(ei,ej,ones(lastindex(ei)),nnodes,nnodes)
    return max.(X,X')
end

#=
requires bins to remain static under edge swaps.
otherwise, silent errors will occur.
two ways to check for this are 
    1. remeasure the desired node characteristic and check that it matches 
    2. check total number of edges in each bin. silent errors can occur if 
        in bins that are incorrectly updated. this will result in later checks
        not checking the correct bins so that the number of edges increases over time. 
=#
"""
    _rewire_binned_edges_degree(E_dict::Dict{Tuple}{Set},nsteps::Int = 100)

given binned edges, perform degree preserving edge swaps within bins. 

This does not check that node_membership is invariant under edge swaps.
Silent errors can occur if this is violated.
"""
function _rewire_binned_edges_degree(
            E_dict::Dict{Tuple}{Set},
            nsteps::Int = 100;
            test::Bool=false
)
    ## setting things up 
    #precompute sampling distribution for bins 
    static_keys = collect(keys(E_dict))

    wv = [length(E_dict[x]) for x in static_keys]
    wv = wv/sum(wv)
    wv = Categorical(wv)
    if test #maintain edges and diagnose errors
        E = Set{Tuple}()
        for key in static_keys
            for edge in E_dict[key]
                push!(E,edge)
            end
        end
        nedges = length(E)
    end
    
    #main loop for edge swaps. 
    #iterate over bins, performing expected number of swaps on each bin 
    failed_swaps = 0
    @showprogress for i=1:nsteps
        bin_ind = rand(wv)
        bin_key = static_keys[bin_ind] #(core src, core dst)    
        curr_bin = E_dict[bin_key]
        (u,v),(s,t) = rand(curr_bin,2) # u,s have same value of node_membership. similarly with v,t

        count=0
        while count<10 && (length(Set([u,v,s,t]))<4 || (u,t) in curr_bin || (s,v) in curr_bin) 
            bin_ind = rand(wv)
            bin_key = static_keys[bin_ind] 
            curr_bin = E_dict[bin_key]
            (u,v),(s,t) = rand(curr_bin,2)    
            count+=1
        end
        if count==10 #skip updates 
            failed_swaps+=1
            continue
        end
        #update bins 
        delete!(curr_bin,(u,v))
        delete!(curr_bin,(s,t))
        delete!(E_dict[reverse(bin_key)],(v,u))
        delete!(E_dict[reverse(bin_key)],(t,s))

        push!(curr_bin,(u,t))
        push!(curr_bin,(s,v))
        push!(E_dict[reverse(bin_key)],(t,u))
        push!(E_dict[reverse(bin_key)],(v,s))

        if test #perform same operations on E
            delete!(E,(u,v))
            delete!(E,(v,u))
            delete!(E,(s,t))
            delete!(E,(t,s))

            push!(E,(u,t))
            push!(E,(t,u))
            push!(E,(s,v))
            push!(E,(v,s))
            if length(E)!=nedges
                return E_dict,E,(u,v),(s,t)
            end
        end 
    end
    @show failed_swaps
    @show nsteps
    return E_dict
end

function _rewire_edges_degree(
            E::Set,
            nsteps::Int = 100
)
    #main loop for edge swaps 
    failed_swaps = 0
    for i=1:nsteps
        (u,v),(s,t) = rand(E,2)

        count=0
        while count<10 && (length(Set([u,v,s,t]))<4 || (u,t) in E || (s,v) in E) 
            (u,v),(s,t) = rand(E,2)
            count+=1
        end
        if count==10 #skip updates 
            failed_swaps+=1
            continue
        end
        #update bins 
        delete!(E,(u,v))
        delete!(E,(s,t))
        
        push!(E,(u,t))
        push!(E,(s,v))
    end
    @show failed_swaps
    @show nsteps
    return E 
end

"""
    _rewire_binned_edges_degree(E_dict::Dict{Tuple}{Set},nsteps::Int = 100)

given binned edges, perform degree preserving edge swaps within bins. 

This does not check that node_membership is invariant under edge swaps.
Silent errors can occur if this is violated.
"""
function _rewire_binned_edges_degree_fast(
            E_dict::Dict{Tuple}{Set},
            nsteps::Int = 100;
            test::Bool=false
)
    ## setting things up 
    #precompute sampling distribution for bins 
    static_keys = collect(keys(E_dict))
    wv = [length(E_dict[x]) for x in static_keys]
    wv = wv/sum(wv)
    wv = Categorical(wv)
    #precompute samples in each bin
    bin_swaps = zeros(Int,lastindex(static_keys))
    for x in rand(wv,nsteps)
        bin_swaps[x]+=1
    end
    #maintain edges to keep nedges fixed
    E = Set{Tuple}()
    for key in static_keys
        for edge in E_dict[key]
            push!(E,edge)
        end
    end
    nedges = length(E)
    
    #main loop for edge swaps. 
    #iterate over bins, performing all swaps  
    failed_swaps = 0
    @showprogress for bin_ind=1:lastindex(static_keys)
        nswaps = bin_swaps[bin_ind]
        bin_key = static_keys[bin_ind] #(core src, core dst)    
        curr_bin = E_dict[bin_key]

        for k=1:nswaps
            (u,v),(s,t) = rand(curr_bin,2) # u,s have same value of node_membership. similarly with v,t

            count=0
            while count<10 && (length(Set([u,v,s,t]))<4 || (u,t) in curr_bin || (s,v) in curr_bin || (u,t) in E || (s,v) in E) 
                (u,v),(s,t) = rand(curr_bin,2)    
                count+=1
            end
            if count==10 #skip updates 
                failed_swaps+=1
                continue
            end
            #update bins 
            #update bins 
            delete!(curr_bin,(u,v))
            delete!(curr_bin,(s,t))
            delete!(E_dict[reverse(bin_key)],(v,u))
            delete!(E_dict[reverse(bin_key)],(t,s))

            push!(curr_bin,(u,t))
            push!(curr_bin,(s,v))
            push!(E_dict[reverse(bin_key)],(t,u))
            push!(E_dict[reverse(bin_key)],(v,s))

            delete!(E,(u,v))
            delete!(E,(v,u))
            delete!(E,(s,t))
            delete!(E,(t,s))

            push!(E,(u,t))
            push!(E,(t,u))
            push!(E,(s,v))
            push!(E,(v,s))
            if length(E)!=nedges
                @show nedges
                @show length(E)
                @show k
                @show bin_ind
                @show count
                return E_dict,E,(u,v),(s,t)
                # @show (u,v,)
            end
        end #end bin swapping in single bin
        #update E_dict
        E_dict[bin_key] = curr_bin
    end
    @show failed_swaps
    @show nsteps
    return E_dict
end


#=
########### TESTING ############


pdensity,nnodes,ngroups = _nested_density_matrix(20,7,0.55,8)
groupsize = Int(nnodes/ngroups)

Random.seed!(8)
A = sparse(stochastic_block_model(pdensity,repeat([round(Int,nnodes/ngroups)],ngroups)))
spy(A)

node_grouping = map(x->floor(Int,(x-1)/groupsize)+1,1:lastindex(A,1))

#test edge binning and reconstruction
E_dict = _get_binned_edges(A,node_grouping)
B = edges_to_graph(E_dict)
norm(A-B) == 0.0

### test rewiring function 
E_dict = _rewire_binned_edges_degree(E_dict,100000);
B = edges_to_graph(E_dict)

#check total edges 
nnz(B)==nnz(A)

#check degrees 
deg_A = vec(sum(A,dims=1))
deg_B = vec(sum(B,dims=1))
norm(deg_A-deg_B)==0

#check sparsity pattern 
spy(A)
spy(B)

#check number of different elements
nnz((A-B).!=0)
nnz(A)
=#