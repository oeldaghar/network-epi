module DiffusionAlgorithms
using SparseArrays
using ProgressMeter
using Printf
using Statistics
using Random
#using PmapProgressMeter
using DataFrames
using DataStructures

import LinearAlgebra.checksquare

#import Compat.LinAlg.checksquare # updates for v0.5

function normout!(A::SparseMatrixCSC{Float64,Int64})
    d = sum(A,dims=2) # sum over rows
    # use some internal julia magic
    for i=1:length(A.nzval)
        A.nzval[i] = A.nzval[i] / d[A.rowval[i]]
    end
    return A
end

"""
- `maxresidvol::Int` - the maximum residual volume considered, if this is negative,
then we treat it as infinite.

Returns
-------
(x::Dict{Int,Float64},r::Dict{Int,Float64},flag::Int)
"""
function weighted_ppr_push(A::SparseMatrixCSC{T,Int}, seed::Int,
    alpha::Float64, eps::Float64, maxpush::Int, dvec::Vector{Int}, maxresidvol::Int) where T

    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval

    n = size(A,1)

    x = Dict{Int,Float64}()     # Store x, r as dictionaries
    r = Dict{Int,Float64}()     # initialize residual
    #Q = Int[]        # initialize queue
    Q = Queue{Int}()
    npush = 0.

    if maxresidvol <= 0
        maxresidvol = typemax(Int)
    end

    rvol = 0

    # TODO handle a generic seed
    r[seed] = 1.
    enqueue!(Q,seed)

    pushcount = 0
    pushvol = 0

    @inbounds while length(Q) > 0 && pushcount <= maxpush
        pushcount += 1
        u = dequeue!(Q)

        du = dvec[u] # get the degree

        pushval = r[u] - 0.5*eps*du
        x[u] = get(x,u,0.0) + (1-alpha)*pushval
        r[u] = 0.5*eps*du

        pushval = pushval*alpha

        for nzi in colptr[u]:(colptr[u+1] - 1)
            pushvol += 1
            v = rowval[nzi]
            dv = dvec[v] # degree of v

            rvold = get(r,v,0.)
            if rvold == 0.
                rvol += dv
            end
            rvnew = rvold + pushval*nzval[nzi]/du

            r[v] = rvnew
            if rvnew > eps*dv && rvold <= eps*dv
                #push!(Q,v)
                enqueue!(Q,v)
            end
        end

        if rvol >= maxresidvol
            return x, r, -2
        end
    end

    if pushcount > maxpush
        return x, r, -1, pushcount
    else
        return x, r, 0, pushcount
    end
end

function weighted_ppr_push_solution(A::SparseMatrixCSC{T,Int}, alpha::Float64, seed::Int, eps::Float64) where T
    maxpush = round(Int,max(1.0/(eps*(1.0-alpha)), 2.0*10^9))
    dvec = sum(A,dims=2)
    return weighted_ppr_push(A,seed,alpha,eps,maxpush,vec(dvec),0)[1]
end

function weighted_local_sweep_cut(A::SparseMatrixCSC{T,Int}, x::Dict{Int,V}, dvec::Vector{Int}, Gvol::Int) where {T,V}
    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval

    n = size(A,1)

    sx = sort(collect(x), by=x->x[2], rev=true)
    S = Set{Int64}()
    volS = 0.
    cutS = 0.
    bestcond = 1.
    beststats = (1,1,1,Gvol-1)
    bestset = Set{Int64}()
    for p in sx
        if length(S) == n-1
            break
        end
        u = p[1] # get the vertex
        #volS += colptr[u+1] - colptr[u]
        volS += dvec[u]

        for nzi in colptr[u]:(colptr[u+1] - 1)
            v = rowval[nzi]
            ew = nzval[nzi]

            if v in S
                cutS -= ew
            else
                cutS += ew
            end
        end
        push!(S,u)
        if cutS/min(volS,Gvol-volS) <= bestcond
            bestcond = cutS/min(volS,Gvol-volS)
            bestset = Set(S)
            beststats = (cutS,min(volS,Gvol-volS),volS,Gvol-volS)
        end
    end
    return bestset, bestcond, beststats
end

function weighted_degree_normalized_sweep_cut!(A::SparseMatrixCSC{T,Int}, x::Dict{Int,V}, dvec::Array{Int}, Gvol::Int) where {T,V}
    colptr = A.colptr
    rowval = A.rowval

    for u in keys(x)
        x[u] = x[u]/dvec[u]
    end

    return weighted_local_sweep_cut(A,x,dvec,Gvol)
end

function weighted_ppr_grow_one(A::SparseMatrixCSC{T,Int},
        seed::Int, alpha::Float64, eps::Float64) where T
    maxpush = round(Int,max(1.0/(eps*(1.0-alpha)), 2.0*10^9))
    dvec = vec(sum(A,dims=2))
    @assert eltype(dvec)==Int
    Gvol = sum(dvec)
    ppr = weighted_ppr_push(A,seed,alpha,eps,maxpush,dvec, 0)[1]
    return weighted_degree_normalized_sweep_cut!(A,ppr,dvec,Gvol)
end

struct setstats
    size::Int
    volume_seed::Int #
    cut::Int
    seed::Int
    volume_other::Int # the volume on the other side
    cond::Float64
    support::Int # total vertices evaluated
    work::Int # total number of edges "pushed"
end

function setstats()
    return setstats(0, 0, 0, 0, 0, 1.0, 0, 0)
end


"""
`weighted_ncp`
--------------

Compute the NCP of a weighted graph. This runs a standard
personalized PageRank-based NCP except on a graph where
the edges are weighted. The weighted edges must be integers
at the moment.

Usage
-----
weighted_ncp(A)
weighted_ncp(A; ...)
- `alpha::Float64` default value 0.99, the pagerank value of alpha
- `maxcond::Float64` ignore clusters with conductance larger than maxcond
- `minsize::Int` the minimum size cluster to consider
- `maxvol::Float64' the maximum volume (if maxvol <= 1, then it's a ratio, if
    maxvol > 1, then it's a total number of edges)
- `maxvisits::Int` don't start a new cluster from a node if it's
  already been in k clusters
Returns
-------
The function returns an Array of setstats
"""

function parallel_weighted_ncp(A::SparseMatrixCSC{Int,Int};
        alpha::Float64=0.99,
        maxcond::Float64=1.0,
        minsize::Int=5,
        maxvol::Float64=1.0,
        maxvisits::Int=10)

    #n = size(A,1)
    n = checksquare(A)
    visits = SharedArray(Int,(n,))

    for i in 1:n
        visits[i] = 0
    end

    eps=1.e-8

    l = ReentrantLock()

    dvec = vec(sum(A,dims=2)) # weighted degree vectors
    @assert eltype(dvec)==Int
    Gvol = sum(dvec)

    ncpdata = pmap(v ->
        begin
            # check if the node has already been handled
            if visits[v] >= maxvisits
                return setstats()
            end

            seed = v

            maxpush = round(Int,max(1.0/(eps*(1.0-alpha)), 2.0*10^9))
            ppr = weighted_ppr_push(A,seed,alpha,eps,maxpush,dvec, 0)[1]
            bestset,bestcond,setnums =
                weighted_degree_normalized_sweep_cut!(A,ppr,dvec,Gvol)

            volseed = setnums[2]
            volother = Gvol-volseed
            if seed ∉ bestset
                volother = volseed
                volseed = Gvol - volseed
            end

            lock(l)
            for u in bestset
                visits[u] += 1
            end
            unlock(l)

            # todo add support/work
            return setstats(length(bestset), volseed, setnums[1], seed, volother, bestcond, 0,0)
        end,
        Progress(n), 1:n)

    return ncpdata
end

function create_ncpdata()
    return DataFrame(seed = Int64[], eps = Float64[], size=Int64[],
                    cond = Float64[], cut = Float64[], volume_seed = Int64[],
                    volume_other = Int64[], ComputeTime = Float64[])
end

function ncp_entry(n, seed, eps, bestset, bestcond, setnums, dt)
    voltrue = setnums[2]
    volseed = setnums[3]
    volother = setnums[4]
    setsize = length(bestset)
    if volother < volseed
        setsize = n - setsize # we are in the complement set
    end

    return [seed, eps, setsize,
            bestcond, setnums[1],
            volseed, volother, dt]
end


function serial_weighted_ncp(A::SparseMatrixCSC{Int,Int};
        alpha::Float64=0.99,
        maxcond::Float64=1.0,
        minsize::Int=5,
        maxvol::Float64=1.0,
        maxvisits::Int=10,
        maxlarge::Int=10,
        largethresh::Float64=0.8,
        epsbig::Float64=1.0e-4)

    n = checksquare(A)

    epsseq = [2,5,10,25,50,100]
    epsvals = [0.01;1.0./(100*epsseq); 1.0./(10000*epsseq); 1.0e-7; 1.0e-8]
    # epsbig = 1.0e-4

    dvec = vec(sum(A,dims=2)) # weighted degree vectors
    @assert eltype(dvec)==Int
    Gvol = sum(dvec)
    meand = mean(dvec)

    ncpdata = create_ncpdata()

    lasteps = false
    for eps=epsvals
        # reset visited
        visits = zeros(n)
        maxpush = round(Int,min(1.0/(eps*(1.0-alpha)), 2.0*10^9))
        if maxpush >= 100*Gvol
            # these are getting too big, let's stop here
            lasteps = true
        end
        vset = randperm(n)
        if maxpush*meand^2 > Gvol
            vset = vset[1:min(ceil(Int64,0.1*n),10^5)]
            if maxpush*meand > Gvol
                vset = vset[1:min(ceil(Int64,0.1*length(vset)),10^4)]
            end
        end
        nlarge = 0
        for v=vset # randomize the order
            seed = v
            if visits[v] >= maxvisits
                continue
            end

            dt = @elapsed ppr = weighted_ppr_push(A,seed,alpha,eps,maxpush,dvec,0)[1]
            dt += @elapsed bestset,bestcond,setnums =
                weighted_degree_normalized_sweep_cut!(A,ppr,dvec,Gvol)

            if length(bestset) <= minsize
                continue
            end
            push!(ncpdata,
                ncp_entry(n, seed, eps, bestset, bestcond, setnums, dt))

            if length(bestset) > largethresh*n
                nlarge += 1
                if nlarge >= maxlarge
                    break
                end
            end

            if eps < epsbig
                for u in bestset
                    visits[u] += 1
                end
            end

            # println(@sprintf("%8.3e  %7i  %8i  %8i  %5.3f  %6.2f",
            #     eps, ncpdata[:size][end], ncpdata[:volume_seed][end], ncpdata[:volume_other][end], bestcond, dt))
        end

        if lasteps
            break
        end
    end

    return ncpdata
end

function serial_weighted_ncp(A::SparseMatrixCSC{Int,Int};
        alpha::Float64=0.99,
        maxcond::Float64=1.0,
        minsize::Int=5,
        maxvol::Float64=1.0,
        maxvisits::Int=10,
        maxlarge::Int=10,
        largethresh::Float64=0.8,
        get_sets::Bool=true,
        epsbig::Float64=1.0e-4)

    n = checksquare(A)

    epsseq = [2,5,10,25,50,100]
    epsvals = [0.01;1.0./(100*epsseq); 1.0./(10000*epsseq); 1.0e-7; 1.0e-8]
    # epsbig = 1.0e-4

    dvec = vec(sum(A,dims=2)) # weighted degree vectors
    @assert eltype(dvec)==Int
    Gvol = sum(dvec)
    meand = mean(dvec)

    ncpdata = create_ncpdata()

    lasteps = false
    ### modify to record best set
    bestsets = Vector{Set{Int64}}()
    ### end modification
    for eps=epsvals
        # reset visited
        visits = zeros(n)
        maxpush = round(Int,min(1.0/(eps*(1.0-alpha)), 2.0*10^9))
        if maxpush >= 100*Gvol
            # these are getting too big, let's stop here
            lasteps = true
        end
        vset = randperm(n)
        if maxpush*meand^2 > Gvol
            vset = vset[1:min(ceil(Int64,0.1*n),10^5)]
            if maxpush*meand > Gvol
                vset = vset[1:min(ceil(Int64,0.1*length(vset)),10^4)]
            end
        end
        nlarge = 0
        for v=vset # randomize the order
            seed = v
            if visits[v] >= maxvisits
                continue
            end

            dt = @elapsed ppr = weighted_ppr_push(A,seed,alpha,eps,maxpush,dvec,0)[1]
            dt += @elapsed bestset,bestcond,setnums =
                weighted_degree_normalized_sweep_cut!(A,ppr,dvec,Gvol)

            if length(bestset) <= minsize
                continue
            end
            push!(ncpdata,
                ncp_entry(n, seed, eps, bestset, bestcond, setnums, dt))
            #modded to store bestset
            if get_sets
                push!(bestsets,bestset)
            end
            #end mod

            if length(bestset) > largethresh*n
                nlarge += 1
                if nlarge >= maxlarge
                    break
                end
            end

            if eps < epsbig
                for u in bestset
                    visits[u] += 1
                end
            end

            # println(@sprintf("%8.3e  %7i  %8i  %8i  %5.3f  %6.2f",
            #     eps, ncpdata[:size][end], ncpdata[:volume_seed][end], ncpdata[:volume_other][end], bestcond, dt))
        end

        if lasteps
            break
        end
    end

    return ncpdata,bestsets
end

function local_serial_weighted_ncp(A::SparseMatrixCSC{Int,Int},seed::Int;
        alpha::Float64=0.99,
        maxcond::Float64=1.0,
        minsize::Int=5,
        maxvol::Float64=1.0,
        maxvisits::Int=10,
        maxlarge::Int=10,
        largethresh::Float64=0.8)

    n = checksquare(A)
    @assert(0<seed<=n,"invalid seed node")

    epsseq = [2,5,10,25,50,100]
    epsvals = [0.01;1.0./(100*epsseq); 1.0./(10000*epsseq); 1.0e-7; 1.0e-8]
    epsbig = 1.0e-4

    dvec = vec(sum(A,dims=2)) # weighted degree vectors
    @assert eltype(dvec)==Int
    Gvol = sum(dvec)
    meand = mean(dvec)

    ncpdata = create_ncpdata()

    lasteps = false
    ### modify to record best set
    bestsets = Vector{Set{Int64}}()
    ### end modification
    for eps=epsvals
        # reset visited
        visits = zeros(n)
        maxpush = round(Int,min(1.0/(eps*(1.0-alpha)), 2.0*10^9))
        if maxpush >= 100*Gvol
            # these are getting too big, let's stop here
            lasteps = true
        end
        # vset = randperm(n)
        vset = [seed]
        # if maxpush*meand^2 > Gvol
        #     vset = vset[1:min(ceil(Int64,0.1*n),10^5)]
        #     if maxpush*meand > Gvol
        #         vset = vset[1:min(ceil(Int64,0.1*length(vset)),10^4)]
        #     end
        # end
        nlarge = 0
        for v=vset # randomize the order
            seed = v
            if visits[v] >= maxvisits
                continue
            end

            dt = @elapsed ppr = weighted_ppr_push(A,seed,alpha,eps,maxpush,dvec,0)[1]
            dt += @elapsed bestset,bestcond,setnums =
                weighted_degree_normalized_sweep_cut!(A,ppr,dvec,Gvol)

            if length(bestset) <= minsize
                continue
            end
            push!(ncpdata,
                ncp_entry(n, seed, eps, bestset, bestcond, setnums, dt))
            #modded to store bestset
            push!(bestsets,bestset)
            #end mod

            if length(bestset) > largethresh*n
                nlarge += 1
                if nlarge >= maxlarge
                    break
                end
            end

            if eps < epsbig
                for u in bestset
                    visits[u] += 1
                end
            end

            # println(@sprintf("%8.3e  %7i  %8i  %8i  %5.3f  %6.2f",
            #     eps, ncpdata[:size][end], ncpdata[:volume_seed][end], ncpdata[:volume_other][end], bestcond, dt))
        end

        if lasteps
            break
        end
    end

    return ncpdata,bestsets
end

function test_delocalization(A::SparseMatrixCSC, alpha::Float64, eps::Float64, dvec::Vector{Int}, Gvol::Int)
    maxpush = round(Int,max(1.0/(eps*(1.0-alpha)), 2.0*10^9))
    ntrials = 5
    for i=1:ntrials
        seed = rand(1:size(A,1))
        ppr,r,flag = weighted_ppr_push(A,seed,alpha,eps,maxpush,dvec,round(Int,Gvol*0.8))
        #@show eps, length(ppr), length(r)
        if flag == -2
            return true
        end

        bestset,bestcond,setnums = weighted_degree_normalized_sweep_cut!(A,ppr,dvec,Gvol)
        #@show eps, setnums, length(bestset)
        if setnums[3] > setnums[4]
            return true
        end
    end

    return false
end

function estimate_delocalization(A::SparseMatrixCSC{Int}, alpha::Float64)

    dvec = vec(sum(A,dims=2))
    Gvol = sum(dvec)

    epscur = 1.0/((1-alpha)*Gvol)
    @assert epscur > eps(1.0)

    rval = test_delocalization(A, alpha, epscur, dvec, Gvol)

    if rval
        epsrange = (epscur, 0.1)
        # already delocalized, increase eps
        for i=1:20
            epscur = epscur*2.0
            if test_delocalization(A, alpha, epscur, dvec, Gvol)
                epsrange = (epscur/2., epscur)
                break
            end
        end
    else
        # decrease eps to delocalize
        epsrange = (eps(1.), epscur)
        for i=1:20
            epscur /= 2.
            if test_delocalization(A, alpha, epscur, dvec, Gvol)
                epsrange = (epscur, epscur*2.)
                break
            end
        end
    end

    # run 5 steps of bisection
    for i=1:5
        mid = epsrange[1]*0.5 + epsrange[2]*0.5
        if test_delocalization(A,alpha,mid, dvec, Gvol)
            # midpoint is delocalized, use upper points
            epsrange = (mid, epsrange[2])
        else
            # could get bigger...
            epsrange = (epsrange[1],mid)
        end
    end

    return epsrange
end

struct NCPProblem
    """ The matrix to run on. """
    A::SparseMatrixCSC{Int,Int}
    dvec::Vector{Int}
    Gvol::Int

    """ The set of alpha values to use. """
    alphas::Vector{Float64}

    """ The largest conductance to consider reporting. """
    maxcond::Float64
    minsize::Int

    """ The maximum volume to consider in a diffusion. """
    maxvol::Float64

    """ The range of eps values to consider. """
    epsrange::Tuple{Float64,Float64}
end

"""
test
"""
function ncp_problem(A::SparseMatrixCSC{Int,Int})
    epsrange = estimate_delocalization(A,0.99) # epsrange[1] is delocalized, epsrange[2] isn't
    dvec = vec(sum(A,dims=2))
    Gvol = sum(dvec)

    return NCPProblem(A, dvec, Gvol, [0.99], 1., 5, 1., epsrange)
end

randlog(a,b) = exp(rand(Float64)*(log(b)-log(a))+log(a))

function random_sample_ncp(p::NCPProblem,N::Int)
    n = size(p.A,1)
    N = min(N,n)
    verts = randperm(n)
    ncpdata = create_ncpdata()

    for i=1:N
        seed = verts[i]
        epsmax = 1/(p.dvec[seed]+1) # just go beyond trivial work
        eps = randlog(p.epsrange[1], epsmax)
        for alpha in p.alphas
            maxpush = round(Int,max(1.0/(eps*(1.0-alpha)), 2.0*10^9))
            dt = @elapsed ppr = weighted_ppr_push(p.A,seed,alpha,eps,maxpush,p.dvec,0)[1]
            dt += @elapsed bestset,bestcond,setnums =
                weighted_degree_normalized_sweep_cut!(p.A,ppr,p.dvec,p.Gvol)

            push!(ncpdata,
                ncp_entry(n, seed, eps, bestset, bestcond, setnums, dt))

            # println(@sprintf("%8.3e  %7i  %8i  %8i  %5.3f  %6.2f",
            #     eps, ncpdata[:size][end],  ncpdata[:volume_seed][end],  ncpdata[:volume_other][end], bestcond, dt))

        end
    end

    return ncpdata
end



#=




function serial_weighted_ncp(A::SparseMatrixCSC{Int,Int};
        alpha::Float64=0.99,
        maxcond::Float64=1.0,
        minsize::Int=5,
        maxvol::Float64=1.,
        maxvisits::Int=10,
        maxlarge::Int=10,
        largethresh::Float64=0.8)
end

=#

end # end module
