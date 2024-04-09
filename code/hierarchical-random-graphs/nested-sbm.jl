# script to generate a nested SBM. Basically the diagonal blocks are dense and get sparser as we go up

#=
density matrix on depth 2 model hierarchy looks like 

[
    p1 p2 p3 p3; 
    p2 p1 p3 p3; 
    p3 p3 p1 p2; 
    p3 p3 p2 p1; 
]

where p1>p2>p3

for the hierarchy 
            O                       - whole network 
    O               O               - bisection
O       O       O       O           - smallest granularity 
=#

using SparseArrays
using LinearAlgebra
using MatrixNetworks
using Graphs
using Distributions
using Plots #for visualizing density matrices

"""
    _get_group(row,col,decay=0.75)

given row and col, figure out what grouping they are in for nested sbm density matrix
"""
function _get_group(row,col,decay=0.75)
    offset = 0.0
    denom = 1.0
    p = 1.0
    i = 1
    maxiter = 5000

    while floor((col+offset)/denom)!=floor((row+offset)/denom) && i<=maxiter
        offset+=denom
        denom*=2
        p*=decay
        i+=1
    end
    return p
end

"""
    _nested_density_matrix(groupsize=10,nlayers=10,decay=0.9,avgd=20)

TBW
"""
function _nested_density_matrix(groupsize=10,nlayers=10,decay=0.9,avgd=20)

    ngroups = 2^nlayers
    nnodes = ngroups*groupsize

    pdensity = zeros(ngroups,ngroups)

    for row=1:ngroups
        for col=1:ngroups
            p = _get_group(row,col,decay)
            pdensity[row,col]=p
        end
    end
    
    #post processing
    #normalize so that rows sum to avgd
    s = sum(pdensity[1,:])
    pdensity*=(avgd/s)

    #enforce constraint on group size. avgd between groups must be feasible
    if any(pdensity.>groupsize-1)
        @warn("enforcing avgd and groupsize")
        pdensity = min.(pdensity,groupsize-1)
    end
    
    return pdensity,nnodes,ngroups
end

"""
    nested_sbm(groupsize,nlayers,avgd=20,decay_factor=0.2)

generates nested SBM with 2^nlayers groups with groupsize nodes in each group.
decay_factor is used to specify the nesting pattern. 
"""
function nested_sbm(groupsize,nlayers,avgd=20,decay_factor=0.2)
    #generate density matrix 
    pdensity,nnodes,ngroups = _nested_density_matrix(groupsize,nlayers,decay_factor,avgd)
    #generate sbm 
    A = sparse(stochastic_block_model(pdensity,repeat([round(Int,nnodes/ngroups)],ngroups)))
    return A
end


############# TESTING ##############
# using Plots
# decay_matrix,_ = _nested_density_matrix(10,1,0.95)
# heatmap(decay_matrix)

# decay_matrix,_ = _nested_density_matrix(10,2,0.95)
# heatmap(decay_matrix)

# decay_matrix,_ = _nested_density_matrix(10,3,0.9)
# heatmap(decay_matrix)

# decay_matrix,_ = _nested_density_matrix(10,5,0.95)
# heatmap(decay_matrix)

#=
#example usage
groupsize = 50
nlayers = 5
decay_factor=0.5
avgd = 20 #should be LESS THAN groupsize or we truncate to min(groupsize-1,avgd)

#generate density matrix 
pdensity,nnodes,ngroups = _nested_density_matrix(groupsize,nlayers,decay_factor,avgd)

@show nnodes 
@show ngroups

#generate sbm 
A = sparse(stochastic_block_model(pdensity,repeat([round(Int,nnodes/ngroups)],ngroups)))
#check average degree 
nnz(A)/lastindex(A,1)
#sparsity patterns 
spy(A)
=#