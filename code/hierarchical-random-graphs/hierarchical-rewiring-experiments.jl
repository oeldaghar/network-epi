using Plots
using Measures
using Random
using ProgressMeter

mainDir = joinpath("/",split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))]...)
dstDir =  joinpath(mainDir,"input","graphs")

include(joinpath(mainDir,"code","graph-io.jl"))
include(joinpath(mainDir,"code","ncp","ncp-acl.jl"))
include(joinpath(mainDir,"code","hierarchical-random-graphs/nested-sbm.jl"))
# include(joinpath(mainDir,"code","hierarchical-random-graphs/hierarchical-rewiring.jl"))


#first pass 
gname = getgnames("study-20","input/graphs/")[1]

#load graph 
A = loadGraph(gname,"input/graphs/")

#perform ncp computations and retrieve sets 
ncp,sets = make_ncp(A,get_sets=true)


sets


#set info 
#set membership (static), set cut, vol (dynamic info)
#node2set 

#edge info 
#edge - set data
    #get from node membership info 

#get candidate quickly 
#way to get edges in a set 
#binned edges by set keys 

# sets
# set2node[set_id] = [nodes in set set_id]
# node2set[node_id] = [sets containing node_id]

# edge2set[(u,v)] = union(node2set[u],node2set[v]) #all sets touched by edge
# set2edge[set_id] -> set2node[set_id] #all nodes in set_id #get their neighbors -> edges that touch this set 

function conductance_sampling(node,node2set,ncp)
    set_inds = collect(node2set[node])

    cond_vals = 1 .-ncp[set_inds,"cond"]
    cond_vals./=sum(cond_vals)

    return set_inds[rand(Categorical(cond_vals))] #sampled set 
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
    ## static data structures
    # set2node = Vector{Set{Int}}()
    # for set_id=1:lastindex(sets)
    #     push!(set2node, sets[set_id])
    # end

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


new_edges = conductance_rewiring(A,ncp,sets,500)
X = edges2graph(new_edges)

#perform ncp computations and retrieve sets 
# newncp,newsets = make_ncp(X,get_sets=true)

# x = map(x->min(newncp.size[x],size(X,1)-newncp.size[x]),1:size(newncp,1))
# y = newncp.cond

# extrema(x)
# extrema(y)

# f = scatter(x,y,xscale=:log10,yscale=:log10,leg=false,
#     xlabel="Size",
#     ylabel="Conductance",
#     guidefontsize=16,
#     tickfontsize=12,
#     xlims=(0.8,maximum(x)*2),
#     ylims=(minimum(y)/8,1.05),
#     top_margin=4Measures.mm)

f,_ = ppr_sample_figure(A)
# f,_ = epi_sample_figure(A)
f

f,_ = ppr_sample_figure(X)
f
# f,_ = epi_sample_figure(X)
# f

# #testing set2edge
# tmp = collect(set2node[1])

# ei,ej,_ = findnz(A[tmp,tmp])
# Set(zip(tmp[ei],tmp[ej]))==set2edge(1,1,set2node,adjlist)

