##
using MatrixNetworks, SparseArrays, LinearAlgebra, Distributions, DelimitedFiles, ProgressMeter, Random, StatsBase, Distributed
##
#=
ctime represents when a node is contagious while rtime represents when
they will recover
=#
module EventEpidemic
#= The idea in the event epidemic.
- use BucketArrays to schedule visiting nodes and node states to decide what action to take.
- we can schedule: moving a node to and from quarantine, infecting a node (S -> E or I)
- if a node is susceptible, look at E.qtimes and E.itime to figure out if we are quarantining or infecting
- if infecting a node, generate incubation and infection times for each neighbor. Then check if these
- are withing the period. If so, add an event to etimes (these are scheduled infection times).
- if node is being infected and quarantined at the same time, quarantine takes precedence so don't attempt to
- infect neighbors.
=#
using DataStructures, SparseArrays, Distributions, Random, DelimitedFiles, MatrixNetworks, StatsBase
import Base.length, Base.setindex!,
  Base.peek, DataStructures.dequeue!, Base.isempty,
  Base.empty!, Base.iterate
mutable struct BucketList
  buckets::Vector{Vector{Int}}
  head::Int
end
""" Create a set structure """
function BucketList(nvals::Int)
  buckets = Vector{Vector{Int}}(undef, nvals)
  for i=1:nvals
    buckets[i] = Vector{Int}()
  end
  bhead = length(buckets)+1
  return BucketList(buckets, bhead)
end
function iterate(bset::BucketList, state=(1, 0))
  # find the next element given state
  i,offset = state
  offset += 1
  # there is better logic nere
  if offset <= length(bset.buckets[i])
    return bset.buckets[i][offset], (i, offset)
  else
    # find the next valid bucket
    i += 1
    offset = 1
    while  i <= length(bset.buckets) && offset > length(bset.buckets[i])
      i += 1
    end
    if i <= length(bset.buckets)
      return bset.buckets[i][offset], (i, offset)
    else
      return nothing
    end
  end
end
function setindex!(bset::BucketList, v::Int, k::Int)
  push!(bset.buckets[v], k)
  if bset.head > v
    bset.head = v
  end
end
function removeindex!(bset::BucketList,v::Int,t::Int)
  """
    function to remove node v from bucket t and update head
  """
  deleteat!(bset.buckets[t],_findfirstindex(bset.buckets[t],v))

  #updating head
  for b=bset.head:length(bset.buckets)
    if length(bset.buckets[b]) > 0
      bset.head = b
      return bset.head
    end
  end
  bset.head = length(bset.buckets)+1
  return bset.head
end
function length(bset::BucketList)
  l = 0
  for b=bset.head:length(bset.buckets)
    l += length(bset.buckets[b])
  end
  return l
end
function peek(bset::BucketList)
  return (bset.buckets[bset.head][1], bset.head)
end
function isempty(bset::BucketList)
  return bset.head > length(bset.buckets)
end
function empty!(bset::BucketList)
  for i=1:length(bset.buckets)
    empty!(bset.buckets[i])
  end
  bset.head = length(bset.buckets)+1
end
function dequeue!(bset::BucketList)
  rval = bset.buckets[bset.head][1]
  popfirst!(bset.buckets[bset.head])
  for b=bset.head:length(bset.buckets)
    if length(bset.buckets[b]) > 0
      bset.head = b
      return rval
    end
  end
  bset.head = length(bset.buckets)+1
  return rval
end
mutable struct SIRData
  beta::Float64 #uniform infection probability
  gamma::Float64 #uniform recovery probability
  log1beta::Float64
  log1gamma::Float64
  inNeighborsA::Vector{Vector{Tuple{Int,Int}}} #adjlist representation of in-neighbors of A: (inneighbor, weight)
  outNeighborsA::Vector{Vector{Tuple{Int,Int}}} #adjlist representation of out-neighbors of A: (outneighbor, weight)
  itime::Vector{Int} # infection time
  ctime::Vector{Int} # contagion time
  rtime::Vector{Int} # recovery time
  snodes::Vector{Bool} # still suceptible nodes
  etimes::BucketList #
  qcap::Int #cap on number of nodes that can be quarantined at once
  qcurr::Vector{Int} #curr num of quarantines
  qnodes::Vector{Bool} # node currently being quarantined
  qtimes::Vector{Int} #time a node last entered quarantine
  ihist::Vector{Vector{Int}} #edge infection times. ihist[v] = incoming infection times for node v
  # outihist::Vector{Vector{Int}}
  # qhist::Vector{Bool}
end
mutable struct SEIRData
  beta::Float64
  gamma::Float64
  a::Float64
  log1beta::Float64
  log1gamma::Float64
  inNeighborsA::Vector{Vector{Tuple{Int,Int}}} #adjlist representation of A
  outNeighborsA::Vector{Vector{Tuple{Int,Int}}} #(v,wᵥ)
  itime::Vector{Int} # infection time
  ctime::Vector{Int} # contagion time
  rtime::Vector{Int} # recovery time
  snodes::Vector{Bool} # still suceptible nodes
  etimes::BucketList
  qcap::Int #cap on number of nodes that can be quarantined at once
  qcurr::Vector{Int} #curr num of quarantines
  qnodes::Vector{Bool}
  qtimes::Vector{Int}
  ihist::Vector{Vector{Int}}
  # qhist::Vector{Bool}
end
function _clear!(E::Union{SIRData,SEIRData})
  """
    resetting the data structure E for a new diffusion
  """
  tmax = typemax(eltype(E.itime))
  fill!(E.itime, tmax)
  fill!(E.ctime, tmax)
  fill!(E.rtime, tmax)
  fill!(E.snodes, true)
  empty!(E.etimes)
  fill!(E.qcurr,zero(Int))
  fill!(E.qnodes, false)
  fill!(E.qtimes, tmax)
  fill!.(E.ihist,tmax) #assumes we already initialized this in the set up
  # fill!(E.qhist, false)
end
function getneighbors(A::SparseMatrixCSC)
  """
      inneighbors adjacency list representation of A
      assumes A[i,j]=1 corresponds to a directed edge from i to j
  """
    inNeighbors = Vector{Vector{Tuple{Int,Int}}}(undef,max(size(A)...))
    ei,ej,ev = findnz(A)
    # rowids = rowvals(A)
    # ev = last(findnz(A))
    @inbounds for i = 1:max(size(A)...)
        inNeighbors[i] = collect(zip(ei[nzrange(A,i)],ev[nzrange(A,i)]))
    end
    return inNeighbors
end
function SIRData(A::SparseMatrixCSC;beta::Float64=0.1,gamma::Float64=0.05,tmax::Int=3500,qcap::Int=0)
  """
    take matrix A, and main epidemic parameters and generate data structure E for epidemic
  """
  nnodes = size(A,1)
  itime = Vector{Int}(undef, nnodes)
  etime = similar(itime)
  rtime = similar(itime)
  snodes = Vector{Bool}(undef, nnodes)
  qnodes = Vector{Bool}(undef, nnodes)
  etimes = BucketList(tmax)
  qcurr = Vector{Int}(undef,tmax)
  qtimes = similar(itime)
  inNeighborsA = getneighbors(A)
  outNeighborsA = getneighbors(sparse(A'))
  ihist = map(w->first.(w),inNeighborsA) # this is somehow faster than a deepcopy
  # outihist = getneighbors(sparse(A'))
  # qhist = similar(snodes)

  E = SIRData(beta, gamma, 1/log1p(-beta), 1/log1p(-gamma),
        inNeighborsA, outNeighborsA, itime, etime, rtime, snodes, etimes, qcap, qcurr, qnodes, qtimes, ihist)#,qhist)#, outihist)
  _clear!(E)
  return E
end
function SEIRData(A::SparseMatrixCSC;beta::Float64=0.1,gamma::Float64=0.05,a::Float64=5.0,tmax::Int=3500,qcap::Int=0)
  """
    take matrix A, and main epidemic parameters and generate data structure E for epidemic
  """
  nnodes = size(A,1)
  itime = Vector{Int}(undef, nnodes) # infected and contagious time
  etime = similar(itime) # exposed time
  rtime = similar(itime) # recovery time
  snodes = Vector{Bool}(undef, nnodes)
  qnodes = Vector{Bool}(undef, nnodes)
  etimes = BucketList(tmax)
  qcurr = Vector{Int}(undef,tmax)
  qtimes = similar(itime)
  inNeighborsA = getneighbors(A)
  outNeighborsA = getneighbors(sparse(A'))
  ihist = map(w->first.(w),inNeighborsA) # this is somehow faster than a deepcopy
  # qhist = similar(snodes)

  E = SEIRData(beta, gamma, a, 1/log1p(-beta), 1/log1p(-gamma),
      inNeighborsA, outNeighborsA, itime, etime, rtime, snodes, etimes, qcap, qcurr, qnodes, qtimes, ihist)#, qhist)
  _clear!(E)
  return E
end
function _findfirstindex(input::Vector{Tuple{Int,Int}},val::Int) #faster than findfirst for our purposes
  """
    find first index i such that input[i]==val.
  """
  inputlen = length(input)
  if inputlen==1 && input[1][1]==val
    return 1
  elseif inputlen > 1
    @inbounds for (i,(v,_)) in enumerate(input)
      if v == val
        return i
      end
    end
    error("val must be in input")
  else
    error("input must be of length>0")
  end
end
function _findfirstindex(input::Vector{Int},val::Int) #faster than findfirst for our purposes
  """
    find first index i such that input[i]==val.
  """
  inputlen = length(input)
  if inputlen==1 && input[1]==val
    return 1
  elseif inputlen > 1
    @inbounds for (i,v) in enumerate(input)
      if v == val
        return i
      end
    end
    error("val must be in input")
  else
    error("input must be of length>0")
  end
end
function _quarantine_node(E::Union{SIRData,SEIRData},t::Int,n::Int,qtime::Int,tmax::Int,neighbors::Bool=true;tdelay::Int=1)
  """
    function to quarantine node n for times t:t+qtime-1
  """
  @inbounds if !E.qnodes[n] && E.qcurr[t]<E.qcap
    E.qnodes[n] = true
    # E.qhist[n] = true
    E.qcurr[t:t+qtime-1] .+= 1
    if E.snodes[n]#schedule Q_S->S
      E.etimes[n] = t+qtime
    else # E,I -> Q -> R so need to update rtime accordingly
      E.rtime[n] = t+qtime
    end
    _adjust_itimes(E,t+1,n,qtime,tmax) #t+1 since this is used after processing scheduled events
  else #adjust qtimes
    E.qtimes[n] = typemax(Int)
  end
  if neighbors #we can still try to quarantine neighbors even if n is already quarantined or if there isnt space now
    if tdelay==0
      for (v,w) in E.outNeighborsA[n] # quarantine
        E.qtimes[v] = t
        _quarantine_node(E,t,v,qtime,tmax,false)
      end
    else
      @inbounds for (v,w) in E.outNeighborsA[n] # schedule quarantine
        if !E.qnodes[v] && E.rtime[v]>t && E.qtimes[v]>t+tdelay #not in Q, not recovered, and v->Q not scheduled
          E.etimes[v] = t+tdelay
          E.qtimes[v] = t+tdelay
        end
      end
    end #end tdelay if-statement
  end #end neighbors if-statement
end
function _myfindmin(input::Vector{Int},excluded_index::Int)
  """
    helper function for _adjust_itimes.
    we need to check that the min doesn't happen at excluded_index.
    if it does, need to return next smallest element
    if there is only one item then return -1,-1 so we know to deschedule a node
  """
  l = length(input)
  if l<2
    return -1,-1
  end
  if excluded_index==1
    m,ind = input[end],l
  else
    m,ind = input[1],1
  end
  @inbounds for (i,v) in enumerate(input)
    if v<m
      if i!=excluded_index
        m,ind = v,i
      end
    end
  end
  return m,ind
end
function _adjust_itimes(E::Union{SIRData,SEIRData},t::Int,n::Int,qtime::Int,tmax::Int)
  """
    adjust the infection times of all nodes that were infected by node n. for each
    of those nodes, we need to find the new node that will pass infection and update
    the scheduling.
  """
  if !E.snodes[n] #already scheduled stuff so need to adjust itimes of outneighbors
    @inbounds for (v,w) in E.outNeighborsA[n]
      #for each node, find inf time by n and compare to min inf time
      if E.snodes[v] # node in S #not necessary but might speed things up?
        itimev = E.itime[v]
        if t<itimev && itimev<=tmax #something was already scheduled
          nind = _findfirstindex(E.inNeighborsA[v],n)
          n_to_v_itime = E.ihist[v][nind]
          if itimev==n_to_v_itime #n could be the only one to infect v
            #find second smallest inf time
            m,ind = _myfindmin(E.ihist[v],nind)
            if m==-1 || m>=tmax #remove event from v
              #remove v from E.etimes.buckets[itimev]
              E.itime[v] = typemax(Int)
              E.ihist[v][nind] = typemax(Int)
              removeindex!(E.etimes,v,itimev)
              # deleteat!(E.etimes.buckets[itimev],_findfirstindex(E.etimes.buckets[itimev],v)) #then update head
            elseif m>itimev #reschedule event for v
              E.itime[v] = m
              E.ihist[v][nind] = m
              removeindex!(E.etimes,v,itimev)
              E.etimes[v] = m
            else #m==itimev so no need to reschedule but need to update ihist
              E.ihist[v][nind] = typemax(Int)
            end
          end
        end
      end
    end #end loop over neighbors
  end
end
function _infect_node(E::SIRData, t::Int, n::Int)
  """
    infect node n at time t. schedule the recovery time and attempt to schedule
    infection times for the neighbors of n
  """
  if E.snodes[n]
    E.snodes[n] = false
    etime = t # exposure time is now for node n
    ctime = etime+1 # forces nodes to be in E for 1 day
    rtime = ctime+(-randexp()*E.log1gamma) #rand(Geometric(E.gamma))
    rtime = ceil(Int,rtime)
    E.itime[n] = t
    E.ctime[n] = ctime
    E.rtime[n] = rtime
    # potentially infect all neighbors
    if !E.qnodes[n] #not quarantined
      @inbounds for (j,(v,w)) in enumerate(E.outNeighborsA[n])
        if E.snodes[v]# not yet handled
          #dt = rand(Geometric(E.beta)) # number of interactions needed for an infection
          # rand(rng::AbstractRNG, d::Geometric) = floor(Int,-randexp(rng) / log1p(-d.p))
          # p = clam(p,w*p)
          pedge = clamp(w*E.beta,1e-6,1-1e-6)
          dt = -randexp() * (1/log1p(-pedge)) #handling edge weight
          if dt <= rtime+1-ctime # then we may infect
            dti = floor(Int,dt)
            itimev = clamp(1, ctime+dti, length(E.etimes.buckets))#forces nodes to be infected for at least a day
            #record inf hist even if not min time inf
            # nind = _findfirstindex(E.inNeighborsA[v],n)
            E.ihist[v][_findfirstindex(E.inNeighborsA[v],n)] = itimev
            if ctime + dti < E.itime[v] # then we would win infection...
              E.itime[v] = itimev
              E.etimes[v] = itimev # set the next event time for this node.
            end
          end
        end
      end #end neighbor loop
    end #end qtime check
  end
end
function _infect_node(E::SEIRData, t::Int, n::Int)
  """
    infect node n at time t. schedule the recovery time and attempt to schedule
    infection times for the neighbors of n
  """
  if E.snodes[n] || E.qtimes[n]==t
    E.snodes[n] = false
    etime = t # exposure time is now for node n
    ctime = etime+ceil(Int, rand(Exponential(E.a)))#forces nodes to be in E for at least 1 day
    rtime = ctime+(-randexp() * E.log1gamma) #rand(Geometric(E.gamma))
    rtime = ceil(Int,rtime)
    E.itime[n] = t
    E.ctime[n] = ctime
    E.rtime[n] = rtime
    # potentially infect all neighbors
    if !E.qnodes[n] #not quarantined
      @inbounds for (j,(v,w)) in enumerate(E.outNeighborsA[n])
        if E.snodes[v] # not yet handled
          #dt = rand(Geometric(E.beta)) # number of interactions needed for an infection
          # rand(rng::AbstractRNG, d::Geometric) = floor(Int,-randexp(rng) / log1p(-d.p))
          pedge = clamp(w*E.beta,1e-6,1-1e-6)
          dt = -randexp() * (1/log1p(-pedge))
          if dt <= rtime+1-ctime # then we may infect
            dti = floor(Int,dt)
            itimev = clamp(1, ctime+dti, length(E.etimes.buckets))#forces nodes to be infected for at least a day
            #record inf hist even if not min time inf. need this for Q_S->S portion
            # nind = _findfirstindex(E.inNeighborsA[v],n)
            E.ihist[v][_findfirstindex(E.inNeighborsA[v],n)] = itimev
            if ctime + dti < E.itime[v] # then we would win infection...
              E.itime[v] = itimev
              E.etimes[v] = itimev # set the next event time for this node.
            end
          end
        end
      end #end neighbor loop
    end #end qtime check
  end
end
function _post_quarantine_infection(E::Union{SIRData,SEIRData}, t::Int, n::Int)
  """
    handles the case n:Q_S->S. Check if any in-neighbors of n are infected and if so, do they infect
    node n after it leaves qurantine.
  """
  for (ind,(v,w)) in enumerate(E.inNeighborsA[n]) #nodes that point to n
    if !E.snodes[v] && E.rtime[v]>t#infected
      olditime = E.ihist[n][ind]
      if olditime<typemax(Int) #v able to pass inf to n
        #newitime = time for v to pass inf to n adjusted for the quarantine period.
        #(olditime-E.ctime[v]) = time for v to pass inf to n
        #max(0,t-20-E.itime[v]) = time steps before quarantine that v failed to infect n
        newitime = olditime-E.ctime[v]-max(0,t-20-E.itime[v])+t
        if newitime<E.rtime[v]  #v infectious
          # ntime = E.itime[n]
          if E.itime[n]<t || (E.itime[n]>t && newitime<E.itime[n])
            E.itime[n] = newitime
            E.etimes[n] = newitime
          end
        end
      end
    end #end v currently infected
  end #end inNeighbors loop
end
function epidemic(E::Union{SEIRData,SIRData}, seed::Int;
  tmax::Int = 3500, q::Bool = false, rseed=-1, tdelay=1, kwargs...)
  """
    main epidemic function. At each time step, dequeue from event times and using
    node state, figure out what action is happening for each node.
  """
  #clearing old data and initializing new epidemic
  _clear!(E)
  if !q
    E.qcap = 0
  end
  if rseed==-1
    Random.seed!(abs(rand(Int64)))
  else
    Random.seed!(rseed)
  end
  _infect_node(E, 0, seed)
  t = 1
  lastactive = t
  nnodes = length(E.snodes)
  ninfs = nnodes-sum(E.snodes) #total infs
  @inbounds while t <= tmax
    # if length(Set(E.etimes.buckets[t]))<length(E.etimes.buckets[t])
    #   println("a node appears $(maximum(counts(E.etimes.buckets[t]))) times at time t=$t")
    # end
    while !isempty(E.etimes) && peek(E.etimes)[2] <= t
      n = dequeue!(E.etimes)
      if E.snodes[n] #currently in S
        if E.qnodes[n] && E.qtimes[n]<=t-20 #Q_S->S #E.qtimes condition for ties of Q_S->S and S->E
          E.qnodes[n] = false
          _post_quarantine_infection(E, t, n)
        elseif E.qtimes[n]==t && !(E.qnodes[n]) #n->Q #!(E.qnodes[n]) condtion for ties with n->I and n->Q
          _quarantine_node(E,t,n,20,tmax,false) #scheduled, so no need to quarantine neighbors
        elseif E.itime[n]==t #S->E,I
          _infect_node(E, t, n)
          lastactive = max(lastactive, E.rtime[n])
        end
      else #E,I->Q
        _quarantine_node(E,t,n,20,tmax,false) #scheduled, so no need to Q neighbors
      end
    end #end processing of scheduled things

    if q #quarantine
      if ninfs>100 #inf detectablity threshold in population
        if E.qcurr[t]<E.qcap #space to quarantine
          th = 1/(E.qcurr[t]+1) # more hesitant to fill spaces when space is low
          @inbounds for c = 1:nnodes
            if E.ctime[c]<=t-1 && E.rtime[c]>=t && !E.qnodes[c]
              if rand()>=th #not a very good condition, should somehow depend on c and avgd
                E.qtimes[c] = t
                _quarantine_node(E,t,c,20,tmax,true;tdelay=tdelay)
                th = 1/(E.qcurr[t]+1)
              end
            end
          end
        end
      else
        ninfs = nnodes - sum(E.snodes)
      end
    end
    t+=1
  end
  return lastactive, E
end
function get_infdata(lastactive::Int,E::Union{SIRData,SEIRData})
  """
    returns vectors of net new infs and new infs at each time step
  """
  infs=zeros(Int,lastactive)
  recovs=zeros(Int,lastactive)
  @inbounds for node=1:length(E.ctime)
    ct,rt = E.ctime[node],E.rtime[node]
    if ct<lastactive
      infs[ct]+=1 #time when infectious
    end
    if rt<lastactive
      recovs[rt+1]+=1 #time when recovered
    end
  end
  netinfs = infs.-recovs
  return netinfs,infs
end
function rtrials(k::Int,E::Union{SIRData,SEIRData},seed::Int;
                rseed::Int=-1,
                tmax::Int=3500,
                q::Bool=true,
                nontrivial::Bool=true,
                tdelay::Int=1,
                inflb::Int=1,
                scount::Union{Vector,SubArray} = zeros(Int,length(E.snodes)))
  """
    runs k trials of epidemic on graph A seeded at seed and returns the newinfs and netnewinfs as new,net
  """
  @assert(inflb<length(E.snodes),"inflb is less than number of nodes in graph")
  if rseed == -1
    rseed = Int64(time_ns())
    Random.seed!(rseed)
  else
    Random.seed!(rseed)
  end
  net,new = Vector{Vector{Int}}(),Vector{Vector{Int}}()
  rseeds = k*seed .+ (randperm(k).-1)
  for i=1:k
    #run trials and store
    Random.seed!(rseeds[i])
    l,E = EventEpidemic.epidemic(E,seed,tmax=tmax,q=q,tdelay=tdelay)
    netinfs,newinfs = EventEpidemic.get_infdata(l,E)
    if nontrivial #must infect some other nodes
      while sum(newinfs)<inflb
        # netinfs,newinfs = EventEpidemic.get_infdata(EventEpidemic.epidemic(E,seed,tmax=tmax,q=q,tdelay=tdelay)...)
        l,E = EventEpidemic.epidemic(E,seed,tmax=tmax,q=q,tdelay=tdelay)
        netinfs,newinfs = EventEpidemic.get_infdata(l,E)
      end
    end
    push!(net,netinfs)
    push!(new,newinfs)
    scount .+= E.snodes
  end
  return new,net#,snodes
end
end # end module


## diffusion functions
using Base.Iterators, ParallelDataTransfer

#new functions that assume that E is mutable
function update_E!(E::Union{EventEpidemic.SIRData,EventEpidemic.SEIRData},beta::Float64,gamma::Float64,qcap::Int)
    #update epidemic parameters
    E.beta = beta
    E.gamma = gamma
    E.qcap = qcap
    E.log1beta = 1/log1p(-beta)
    E.log1gamma = 1/log1p(-gamma)
end

update_E!(E::Union{EventEpidemic.SIRData,EventEpidemic.SEIRData},beta::Float64) = update_E!(E,beta,0.05,0)

function epidemic_diffusion(b::Float64,gamma::Float64,qcap::Int,v::Int,ktrials::Int;
                    tmax::Int=3500,tdelay::Int=1,inflb::Int=1) # assumes that E is already loaded in memory
    #update paramaters
    update_E!(E,b,gamma,qcap)
    #pick correct column of snodes to use
    sview = view(scounts,:,qdict[qcap])

    #run diffusion
    if qcap == 0 #should be faster when q=false
      res = EventEpidemic.rtrials(ktrials,E,v,rseed=abs(rand(Int)),tmax=tmax,q=false,tdelay=tdelay,inflb=inflb,scount=sview)
      return res
    else
      res = EventEpidemic.rtrials(ktrials,E,v,rseed=abs(rand(Int)),tmax=tmax,q=true,tdelay=tdelay,inflb=inflb,scount=sview)
      return res
    end
end

function distributed_diffusion(b::Float64,gamma::Float64,qcaps::Vector{Int},vs::Vector{Int},ktrials::Int;
                                a::Float64=5.0,tdelay::Int=1,inflb::Int=1,tmax::Int = 3500,method::String="sir")

    #send A everywhere->build E locally->process parameters
    # @everywhere @eval A=$A

    @everywhere @eval b,gamma,tmax,inflb,tdelay,ktrials,a = $b,$gamma,$tmax,$inflb,$tdelay,$ktrials,$a
    #put array on each thread and just add them up at the end
    nq = length(qcaps)
    @everywhere @eval nq = $nq
    #make dict to make things easier
    qdict = Dict{Int,Int}()
    for (i,val) in enumerate(qcaps)
      qdict[val] = i
    end
    @everywhere @eval qdict = $qdict


    if lowercase(method)=="sir"
      @everywhere @eval E = Main.EventEpidemic.SIRData(A,beta=b,gamma=gamma,tmax=tmax,qcap=0)
    elseif lowercase(method)=="seir"
      @everywhere @eval E = Main.EventEpidemic.SEIRData(A,beta=b,a=a,gamma=gamma,tmax=tmax,qcap=0)
    else
      @error("method must be either sir or seir")
    end
    # println(typeof(E)) #for testing 
    scounts = zeros(Int,size(A,1),length(qcaps))
    @everywhere @eval scounts = $scounts

    #prepping parameters
    # @everywhere tfun((v,qcapacity,beta)) = epidemic_diffusion(beta,gamma,qcapacity,v,ktrials,tdelay=tdelay,inflb=inflb,tmax=tmax)
    ps = vcat(Base.Iterators.product(vs,qcaps,b)...)

    #diffusion
    # for worker in workers()
    #   println("worker $worker: sum is $(sum(fetch(remotecall(x->x,workers()[1],scounts)))))")
    # end
    res = @showprogress pmap(x->epidemic_diffusion(x[3],gamma,x[2],x[1],ktrials,tdelay=tdelay,inflb=inflb,tmax=tmax),ps)
    # println("after pmap call")
    # for worker in workers()
    #   println("worker $worker: sum is $(sum(fetch(remotecall(x->x,workers()[1],scounts)))))")
    # end
    return res
end


# ( base graph triangles * beta ) = (rewired triangles * new beta)
function _find_triangle_beta(num_tris,beta,num_newtris)
  return (num_tris/num_newtris)*beta
end

function qdiffusion1(nnodes::Int,
  ktrials::Int,
  gnames::Vector{String};
  betas::Vector{Float64}=collect(0.025:0.025:0.1),
  gamma::Float64=0.05,
  gpath::String,
  dst::String="diffusion-data/",
  qpercents::Union{Vector{Int},Vector{Float64}}=collect(0:15),
  qsizes::Vector{Int}=Vector{Int}(),
  nbins::Int=10,
  tmax::Int=3500,
  method::String="sir",
  tdelay::Int=1,
  inflb::Int=2,
  nodesampletype::String="degree",
  node_rseed::Int=71)


  sort!(qpercents)
  minq,maxq = extrema(qpercents)
  @assert(0<=minq && maxq<=100)
  # A = loadGraph(gnames[1],gpath,sym=true)
  @everywhere @eval A = loadGraph($(gnames[1]),$gpath,sym=true)
  @everywhere @eval A = max.(A,$commonNeighbors(A))#for triangle rewiring diffusions
  # @everywhere @eval A=$A
  #read graph,symmetrize and kill diagonal
  if length(qsizes)==0
    qsizes = [div(k*size(A,1),100) for k in qpercents]
  end

  @inbounds for (i,g) in enumerate(gnames)
    println("working on graph: $g")
    if g != gnames[1] #dont reload first graph
      @everywhere @eval A = loadGraph($(gnames[1]),$gpath,sym=true)
      @everywhere @eval A = max.(A,$commonNeighbors(A))
    end

    #node sample
    if nodesampletype=="degree"
      ns = stratified_nodesample(vec(sum(Int,A;dims=2)),nnodes,nbins = nbins,rseed=node_rseed)
    elseif nodesampletype=="uniform"
      Random.seed!(node_rseed)
      ns = Distributions.sample(1:size(A,1),nnodes,replace=false)
    elseif nodesampletype=="kcore"
      ns = stratified_nodesample(vec(corenums(A)[1]),nnodes,nbins = nbins,rseed=node_rseed)
    end
    @assert(all(0 .<betas.<1),"must have 0<b<1 for all values of beta")

    ##epidemic simulation
    for (i,b) in enumerate(betas)
      pname="$(betas[i])-$gamma"
      # scounts = zeros(Int,size(A,1),length(qsizes))

      println("$(uppercase(method)) parameter set $i of $(length(betas))")
      # res = distributed_diffusion(A,b,gamma,qsizes,ns,ktrials,tmax=tmax,method=method)
      res = distributed_diffusion(b,gamma,qsizes,ns,ktrials,tmax=tmax,method=method)
      tinfs = vcat(map(x->x[1],res)...) #tinfs
      cinfs = vcat(map(x->x[2],res)...) #cinfs

      #savedata
      open(dst*"tinfs-$(g[1:end-5])-$method-"*pname*".txt","w") do io
        writedlm(io,tinfs,",")
      end
      open(dst*"cinfs-$(g[1:end-5])-$method-"*pname*".txt","w") do io
        writedlm(io,cinfs,",")
      end

      # accumlate scounts on parent
      for worker in workers()
        scounts .+= @getfrom worker scounts
      end #TODO improve to O(log(nworkers)) via summing over pid,2*pid and recursing

      #save scounts
      open(dst*"scounts-$(g[1:end-5])-$method-"*pname*".txt","w") do io
        writedlm(io,scounts,",")
      end
      # CSV.write(dst*"scounts-$(g[1:end-5])-sir-"*pname*".txt",scounts)#alt
    end #end beta loop
  end #end graph loop
end

function triangle_qdiffusion(nnodes::Int,
  ktrials::Int,
  gnames::Vector{String};
  betas::Vector{Float64}=collect(0.025:0.025:0.1),
  gamma::Float64=0.05,
  gpath::String,
  dst::String="diffusion-data/",
  qpercents::Union{Vector{Int},Vector{Float64}}=collect(0:15),
  qsizes::Vector{Int}=Vector{Int}(),
  nbins::Int=10,
  tmax::Int=3500,
  method::String="sir",
  tdelay::Int=1,
  inflb::Int=2,
  nodesampletype::String="degree",
  node_rseed::Int=71)


  sort!(qpercents)
  minq,maxq = extrema(qpercents)
  @assert(0<=minq && maxq<=100)
  # A = loadGraph(gnames[1],gpath,sym=true)
  @everywhere @eval A = loadGraph($(gnames[1]),$gpath,sym=true)
  @everywhere @eval A = max.(A,$commonNeighbors(A))#for triangle rewiring diffusions


  if length(qsizes)==0
    qsizes = [div(k*size(A,1),100) for k in qpercents]
  end

  #cache num of hyperedges for base network
  base_num_tris = sum(A)

  @inbounds for (i,g) in enumerate(gnames)
    println("working on graph: $g")
    if g != gnames[1] #dont reload first graph
      @everywhere @eval A = loadGraph($(gnames[1]),$gpath,sym=true)
      @everywhere @eval A = max.(A,$commonNeighbors(A))
    end

    #node sample
    if nodesampletype=="degree"
      ns = stratified_nodesample(vec(sum(Int,A;dims=2)),nnodes,nbins = nbins,rseed=node_rseed)
    elseif nodesampletype=="uniform"
      Random.seed!(node_rseed)
      ns = Distributions.sample(1:size(A,1),nnodes,replace=false)
    elseif nodesampletype=="kcore"
      ns = stratified_nodesample(vec(corenums(A)[1]),nnodes,nbins = nbins,rseed=node_rseed)
    end
    @assert(all(0 .<betas.<1),"must have 0<b<1 for all values of beta")

    ##epidemic simulation
    for (i,b) in enumerate(betas)
      pname="$b-$gamma"

      #get new beta to keep (total triangles*beta) fixed from original to rewired network
      newb = _find_triangle_beta(base_num_tris,b,sum(A))
      newb = clamp(newb,1e-6,1-1e-6) 

      println("$(uppercase(method)) parameter set $i of $(length(betas))")
      # res = distributed_diffusion(A,b,gamma,qsizes,ns,ktrials,tmax=tmax,method=method)
      res = distributed_diffusion(newb,gamma,qsizes,ns,ktrials,tmax=tmax,method=method,tdelay=tdelay,inflb=inflb)
      tinfs = vcat(map(x->x[1],res)...) #tinfs
      cinfs = vcat(map(x->x[2],res)...) #cinfs

      #savedata
      open(dst*"tinfs-$(g[1:end-5])-$method-"*pname*".txt","w") do io
        writedlm(io,tinfs,",")
      end
      open(dst*"cinfs-$(g[1:end-5])-$method-"*pname*".txt","w") do io
        writedlm(io,cinfs,",")
      end

      # accumlate scounts on parent
      for worker in workers()
        scounts .+= @getfrom worker scounts
      end #TODO improve to O(log(nworkers)) via summing over pid,2*pid and recursing

      #save scounts
      open(dst*"scounts-$(g[1:end-5])-$method-"*pname*".txt","w") do io
        writedlm(io,scounts,",")
      end
      # CSV.write(dst*"scounts-$(g[1:end-5])-sir-"*pname*".txt",scounts)#alt
    end #end beta loop
  end #end graph loop
end

#general functions that help with io, sampling, etc..
function stratified_nodesample(nodemap::Vector,nsamples::Int;nbins::Int=5,rseed::Int=-1)
  """take a stratifed sample of nodes depending on the nodemap (typically degree of corenumber)"""
  if rseed>=0
    Random.seed!(rseed)
  else
    Random.seed!(time_ns())
  end

  bounds = quantile(nodemap,1/nbins:1/nbins:1.0)
  binsamplesize = div(nsamples,nbins)
  offset = nsamples-binsamplesize*nbins

  cmap = countmap(bounds)
  bounds = unique(bounds)
  res=Vector{Int}()

  if length(bounds)<nbins
      p=sortperm(collect(keys(cmap)))
      bweights=collect(values(cmap))[p]
      nbins=length(bounds)
  else
      bweights = ones(nbins)
  end
  #split nodes into bins
  bins = Dict{Int,Vector{Int}}()
  for i = 1:length(bounds)
      bins[i] = Vector{Int}()
  end

  for (i,v) in enumerate(nodemap)
    t=1
    while t<=nbins && v>bounds[t]
      t+=1
    end
    push!(bins[t],i)
  end

  for k=1:nbins
      if offset>0
          samplesize = Int(binsamplesize*bweights[k]+min(bweights[k],offset))
          offset-=min(bweights[k],offset)
      else
          samplesize = Int(binsamplesize*bweights[k])
      end
      push!(res,Distributions.sample(bins[k],replace=false,samplesize)...)
  end
  return res
end

function commonNeighbors(A::SparseMatrixCSC)
  @assert issymmetric(A)
  P = deepcopy(A)
  fill!(nonzeros(P),0)
  rowval = rowvals(A)
  neighs = zeros(Bool, size(A,1))
  for j=1:size(A,1)
    # index neighbors
    for nzj in nzrange(A, j)
      neighs[rowval[nzj]] = true
    end
    # for each edge out of this node...
    for nzj in nzrange(A, j)
      i = rowval[nzj]
      score = 0
      for nzi in nzrange(A, i)
        w = rowval[nzi]
        if neighs[w] == true
          # for each neighbor of i (w) that is
          # also a neighbor of j neighs[w] = true
          # increase the score
          score += 1
        end
      end
      nonzeros(P)[nzj] = score
    end
    # reset indicies
    for nzj in nzrange(A, j)
      neighs[rowval[nzj]] = false
    end
  end
  return P
end

