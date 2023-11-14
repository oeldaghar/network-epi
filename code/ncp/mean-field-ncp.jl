"""
    X is a node-by-timestep matrix with the probability that node i is
    infected at time j given that X was not infected at any previous time-step
    Y is a node-by-timestep matrix with the probability that node i is
    infected at time j
    This function is designed to test the model in the most simple scenario possible

    t is the current timestep, we assume that t > 1 and that all information
    from previous t is specified.

    p is the infection probability
    gamma is the recovery probability
"""
# A = 
using LinearAlgebra
using MatrixNetworks
using SparseArrays

#allows for execution from command line as well as an ide so files can be run modularlly
mainDir = joinpath(split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))])

include(joinpath(mainDir,"code/fast-diffusion.jl"))
include(joinpath(mainDir,"code/graph-io.jl"))

# function evolve_simple!(X::Matrix, Y::Matrix, t::Int, p::Real, gamma::Real, A::AbstractMatrix)
#   @assert(t > 1)
#   @assert(t <= size(X,2) && t <= size(Y,2))
#   logy = log.(1 .- p.*Y[:,t-1])
#   X[:,t] .= (1 .- exp.(A*logy)).*(1 .- vec(sum(@view X[:,1:t-1];dims=2)))# probability that at least one node gave us infection
#   Y[:,t] .= X[:,t]
#   gamma_s = 1.0 # gamma raised to the power s
#   for s=t-1:-1:1
#     gamma_s *= (1-gamma) # not recovered
#     Y[:,t] .+= gamma_s .* X[:,s]
#   end
# end
function evolve_simple!(x::Vector{Float64},
                        y::Vector{Float64},
                        z::Vector{Float64},
                        t::Int,
                        beta::Float64,
                        gamma::Float64,
                        A::SparseMatrixCSC)
  #at initialization, cannot have y==0 or else iterates don't change
  #x = Pr(i gets infected at current time step)
  #y = Pr(i infectious at current time step)
  #z = Pr(i has been infected up to current time step)

  for k = 1:t  
    logy = log.(1 .- beta.*y)
    x .= (1 .- exp.(A*logy)).*(1 .- z)# probability that at least one node gave us infection
    y .= x .+(1-gamma).*y
    z .= z .+ x
  end
  return x,y,z
end


"""
  smooth_epidemic(l,E)
    take a partial (or complete) epidemic and return a vector indicating node ranking
    in particular, recoveries and future infections are ignored. 
""" 
function smooth_epidemic(l::Int,
                          itimes::Vector{Int},
                          timesteps::Int,
                          A::SparseMatrixCSC)
  #zero out future infection history
  x = zeros(Float64,max(size(A)...))
  #make all infected nodes contagious
  y = Float64.(itimes.<=l)
  #all persons infected since t=0
  z = deepcopy(y)

  return evolve_simple!(x,y,z,timesteps,E.beta,0.0,A)
end 

smooth_epidemic(l::Int,E::Union{EventEpidemic.SEIRData,EventEpidemic.SIRData},timesteps::Int,A::SparseMatrixCSC) = smooth_epidemic(l,E.itime, timesteps,A)


#w::Vector{Float64} = normalized node scores (seed node -> 1)
function smooth_epidemic(w::Vector{Float64},beta::Float64,A::SparseMatrixCSC,timesteps::Int=10)
  #zero out future infection history
  x = zeros(Float64,max(size(A)...))
  #use node scores for Pr(node contagious). 
  # y = deepcopy(w)
  #all persons infected since t=0
  # z = deepcopy(w)

  return evolve_simple!(x,w,w,timesteps,beta,0.0,A)
end



#some ideas to try out 
#   smoothing out the epidemic
#     run multiple epidemics and use how many times a node got infected as important (current)
#     at tend, evolve epidemic a few steps combine with E.itime as ranking 
#   looking at difference in predictions 
#     partial epidemic vs predicted epidemic. difference = importance 
#     

using Plots
ENV["GKSwstype"] = "100"

#set things up 
gname = getgnames("lfr","input/graphs/")[end]
A = loadGraph(gname,"input/graphs/")
E = EventEpidemic.SEIRData(A,beta=0.5,tmax=10000);
n = size(A,1)

#smoothing out many epidemics 
seed_node = rand(1:lastindex(E.snodes))
z = zeros(Int,lastindex(E.snodes));
z1 = zeros(Int,lastindex(E.snodes));
z2 = zeros(Float64,lastindex(E.snodes));
curr_failures = 0
kfailures = 6
target_infs = 10000

betas = vcat(1e-2:1e-2:1e-1,2e-1:1e-1:5e-1);
update_E!(E,0.3)

@time for trial_num = 1:30#lastindex(betas)
  #increase beta (infection goes further)
  # update_E!(E,betas[trial_num])
  # println("working on trials $trial_num")
  # l,E = EventEpidemic.partial_epidemic(E,seed_node,target_infs)
  l,E = EventEpidemic.epidemic(E,seed_node)
  netinfs, newinfs = EventEpidemic.get_infdata(l, E)
  # while  sum(newinfs) < target_infs && curr_failures<kfailures
    # curr_failures +=1 
    # #l,E = EventEpidemic.partial_epidemic(E,seed_node,target_infs)
    # l,E = EventEpidemic.epidemic(E,seed_node)
    # netinfs, newinfs = EventEpidemic.get_infdata(l, E)
  # end
  @show sum(newinfs)

  # weight_fcn(x) = E.snodes[x] ? 0.0 : 1/(1+E.itime[x])

  # if curr_failures<kfailures
  z.+=E.snodes  

  z1.+=min.(E.itime,l+1)

  z2 .+= weight_fcn.(1:lastindex(E.snodes))

  
  # end
end

nnz(A)/size(A,1)

curr_failures

z1[seed_node]
z[seed_node]
z2[seed_node]

#weighted by num of times a node was infected 
p = sweepcut(A,lastindex(betas) .- z);
min(argmin(p.conductance),size(A,1)-argmin(p.conductance))
plot(eachindex(p.conductance),p.conductance,
      xscale=:log10,yscale=:log10,ylims=(1e-3,1),leg=false)
extrema(p.conductance)

#weighted by infection time 
p = sweepcut(A,- z1);
min(argmin(p.conductance),n-argmin(p.conductance))
plot(eachindex(p.conductance),p.conductance,
      xscale=:log10,yscale=:log10,ylims=(1e-3,1),leg=false)
extrema(p.conductance)

# weighted by infection time for infected nodes
p = sweepcut(A,z2);
min(argmin(p.conductance),n-argmin(p.conductance))
plot(eachindex(p.conductance),p.conductance,
      xscale=:log10,yscale=:log10,ylims=(1e-3,1),leg=false)
extrema(p.conductance)


n


#smoothing out tail 
w = smooth_epidemic(z2./15,0.01,A,10)[3]
p = sweepcut(A,w);
argmin(p.conductance)
plot(eachindex(p.conductance),p.conductance,
      xscale=:log10,yscale=:log10,ylims=(1e-2,1),leg=false)

# plot(collect(10:10:lastindex(p.conductance)),p.conductance[10:10:lastindex(p.conductance)],
        # xscale=:log10,yscale=:log10,ylims=(1e-2,1),leg=false)
