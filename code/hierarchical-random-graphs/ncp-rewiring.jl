#=
functions for rewiring to approximately preserve the NCP.

note that this method is sensitive to the number of edges rewired. 

=#
using Plots
using Measures
using Random
using ProgressMeter

mainDir = joinpath("/",split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))]...)
dstDir =  joinpath(mainDir,"input","graphs")

include(joinpath(mainDir,"code","graph-io.jl"))
include(joinpath(mainDir,"code","ncp","ncp-acl.jl"))
# include(joinpath(mainDir,"code","hierarchical-random-graphs/nested-sbm.jl"))
# include(joinpath(mainDir,"code","hierarchical-random-graphs/hierarchical-rewiring.jl"))


function compute_set_cond(S::T,A::SparseMatrixCSC) where T
    #
    S_indicator = zeros(Bool,lastindex(A,1))
    for node in S
        S_indicator[node] = true
    end
    
    cut = 0
    vol = 0
    row_ids = rowvals(A)
    for node in S
        for neighbor in row_ids[nzrange(A,node)] 
            vol+=1
            cut+=1-S_indicator[neighbor]
        end
    end
    vol = min(vol,nnz(A)-vol)
    return cut/vol 
end

function conductance_sampling(node,node2set,ncp)
    set_inds = collect(node2set[node])

    #inverts 
    # cond_vals = 1 .-ncp[set_inds,"cond"]
    # cond_vals = cond_vals.^25    
    # cond_vals./=sum(cond_vals)
    # return set_inds[rand(Categorical(cond_vals))] #sampled set 

    return set_inds[argmin(ncp[set_inds,"cond"])]
end
# dynamic data structures
function edge2set(src,dst,node2set)
    return union(node2set[src],node2set[dst])
end

function _get_neighbor_list(A::SparseMatrixCSC)
    neighs = Vector{Set{Int}}(undef,lastindex(A,1))
    rows = rowvals(A)
    for col = 1:lastindex(A,1)
        neighs[col] = Set(rows[nzrange(A,col)])
    end
    return neighs
end

function _get_edgelist(A::SparseMatrixCSC)
    ei,ej,ev = findnz(A)
    return Set(zip(ei,ej))
end

function set2edge(src_set,dst_set,set2node,adjlist)  
    set_edges = Set{Tuple{Int,Int}}()

    #nodes in dst
    dst_nodes = set2node[dst_set]

    #find neighbors of first set that lie in the second set 
    for src_node in set2node[src_set]
        for endpoint in adjlist[src_node]
            if endpoint in dst_nodes
                push!(set_edges,(src_node,endpoint))
            end
        end
    end
    return set_edges
end

function conductance_rewiring(A,ncp,set2node,nsteps=100)
    #cache node2set and set2node memberships
    node2set = Vector{Set{Int}}(undef,lastindex(A,1))
    for node_id=1:lastindex(A,1)
        node2set[node_id] = Set{Int}()
    end
    for (set_id,set) in enumerate(set2node)
        for node in set
            push!(node2set[node],set_id)
        end
    end

    # dynamic data structures 
    adjlist = _get_neighbor_list(A)
    edges = _get_edgelist(A)

    ## main loop 
    failed_swaps = 0
    @showprogress for k=1:nsteps 
        #sample edge 
        (u,v) = rand(edges)

        #sample a set for each endpoint
        u_set = conductance_sampling(u,node2set,ncp)
        v_set = conductance_sampling(v,node2set,ncp)

        #find similar edges 
        tmp_edges = set2edge(u_set,v_set,set2node,adjlist)

        #rewire
        (s,t) = rand(tmp_edges)
        attempt = 0
        while attempt<10 && (length(Set([u,v,s,t]))<4 || (u,t) in edges || (v,s) in edges)
            (s,t) = rand(tmp_edges)
            attempt+=1
        end

        if attempt==10 #skip updates 
            failed_swaps+=1
            continue
        end
        #update edges 
        delete!(edges,(u,v))
        delete!(edges,(v,u))
        delete!(edges,(s,t))
        delete!(edges,(t,s))

        push!(edges,(u,t))
        push!(edges,(t,u))
        push!(edges,(s,v))
        push!(edges,(v,s))

        #update adjlist
        delete!(adjlist[u],v)
        delete!(adjlist[v],u)
        delete!(adjlist[s],t)
        delete!(adjlist[t],s)

        push!(adjlist[u],t)
        push!(adjlist[v],s)
        push!(adjlist[s],v)
        push!(adjlist[t],u)
    end
    return edges 
end

function edges2graph(edges)
    ei,ej = Vector{Int}(),Vector{Int}()
    for (u,v) in edges
        push!(ei,u)
        push!(ej,v)
    end
    nnodes = max(maximum(ei),maximum(ej))
    X = sparse(ei,ej,ones(lastindex(ei)),nnodes,nnodes)
    return max.(X,X')
end

function preprocess_sets!(ncp,sets::Vector{Set{T}},A::SparseMatrixCSC,max_set_size::Int=250) where T
    #filter sets and ncp 
    set_inds = ncp.size.<=max_set_size
    s = Set{Int}()
    for (i,x) in enumerate(sets)
        if set_inds[i]
            union!(s,x)
        end
    end
    ncp = ncp[set_inds,:]
    sets = sets[set_inds]

    #add a new set containing missed nodes and update data structures 
    if length(s)<lastindex(A,1)
        other_nodes = setdiff(1:lastindex(A,1),s)
        set_size = Float64(length(other_nodes))
        set_cond = compute_set_cond(s,A)
        push!(ncp,[-1.0,-1.0,set_size,set_cond,-1.0,-1.0,-1.0,-1.0])
        push!(sets,Set(other_nodes))
    end
    return ncp,sets
end
#=
##### Testing ######


#use small graph for quick testing 
gname = getgnames("geometric","input/graphs/")[1]
g = gname[1:end-5]

# h = "ncp-rewired-100.0-$gname"

fpath = joinpath(mainDir,"input/graphs/",gname)
# MatrixNetworks.readSMAT(fpath)

#load graph 
A = loadGraph(gname,"input/graphs/")
X = loadGraph(h,"pipeline/graphs/")

# #load precomputed NCP and sets 
ncploc = joinpath(mainDir,"pipeline/data/$g/ncpdata/")
ncp,_,sets = load_ncpdata(gname,ncploc)
length(sets)

##alternatively, do ncp computations
# ncp,sets = make_ncp(A,get_sets=true)

histogram(length.(sets))
#preprocess sets 
ncp,sets = preprocess_sets!(ncp,sets,A,5000)
length(sets)

# @showprogress for x in sets
#     inds = collect(x)
#     if size(largest_component(A[inds,inds])[1],1)!=lastindex(inds)
#         @show x 
#     end
# end


# using Profile
new_edges = conductance_rewiring(A,ncp,sets,50000)
X = edges2graph(new_edges)
nnz(A)==nnz(X)
#perform ncp computations and retrieve sets 

f,_ = ppr_sample_figure(A)
# f,_ = epi_sample_figure(A)
f


f,_ = ppr_sample_figure(X)
# f,_ = epi_sample_figure(A)
f


=#


# function ncp_rewire(A,sets,nsteps)
#     #preprocess sets. 
#     filter!(x->set_size<=250,sets)
#     #put all nodes that don't have a set into a single large set
#     new_set = setdiff(1:nnodes,vcat(sets.nodes...))
#     push!(sets,new_set)
#     #end preprocessing

#     for i=1:nsteps
#         #sample a random edge
#         u,v = rand(edges)
#         #sample/lookup set for each endpoint of edge 
#         #currently, this is the minium conductance set containing that node
#         uset,vset = setlookup(u),setlookup(v)
#         #find similar edges that cross those sets in the same manner
#         similar_edges = get_edges(uset,vset)
#         #sample a candidate for swapping
#         s,t = rand(similar_edges)
#         #perfrom edge swap 
#         swap_edges((u,v),(s,t))
#         #update data structures for edges and adjlist
#         update!(edges,(u,v),(s,t))
#         update!(adjlist,(u,v),(s,t))
#     end
#     return edges
# end


#=
f,_ = ppr_sample_figure(A,500)
f


#for each node, find size of min conductance set containing it 
node2set = Vector{Set{Int}}(undef,lastindex(A,1))
for node_id=1:lastindex(A,1)
    node2set[node_id] = Set{Int}()
end
for (set_id,set) in enumerate(sets)
    for node in set
        push!(node2set[node],set_id)
    end
end

# dynamic data structures 
adjlist = _get_neighbor_list(A)
edges = _get_edgelist(A)

cached_set2edge = Dict{Tuple{Int,Int},Set{Int}}

ei,ej,_ = findnz(A)

s = Vector{Int}() 
@showprogress for i=1:lastindex(ei)
    push!(s,length(edge2set(ei[i],ej[i],node2set)))
end

nnz(A)
sum(s)/size(A,1)

extrema(length.(node2set))


wv = 1 .-ncp.cond

wv./=sum(wv)
cond_dist = Categorical(wv)


Distributions.sample(1:lastindex(wv),wv)


sets[inds]

histogram(ncp.cond[inds])


# sample(1:lastindex(wv),weights(wv),500,replace=false)

=#