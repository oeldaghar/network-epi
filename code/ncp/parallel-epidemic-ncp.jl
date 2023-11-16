#parallel_epidemic_ncp
#initialize graph and epidemic data structure
#pmap(beta,trial) -> run diffusion and return relevant results
#assumes A and E are loaded locally on processor
using Distributed, ProgressMeter

using MatrixNetworks, LinearAlgebra, SparseArrays, DataFrames

using DelimitedFiles
using CSV 


function _log_rand_unif(lower,upper)
  
  #transform to log scale 
  log_lower = log10(lower)
  log_upper = log10(upper)

  #sample uniformly in [log_lower,log_upper]
  pts = (log_upper-log_lower).*rand(2).+log_lower

  #transform back to original scale and make integers
  pts = round.(Int,(10.0).^pts)
  return pts
end

function _subsample2(data,splits=[1;100;1000;10000;100000])
  maxlen = lastindex(data)-1

  for i=1:lastindex(splits)
    if splits[i]>lastindex(data)-1
      splits = splits[1:i-1]
      break
    end
  end
  
  inds = Vector{Int}()
  # bounds = []
  for i=2:lastindex(splits)
    #subsample bounds 
    # lowerbound,upperbound = extrema(rand(splits[i-1]:splits[i],2))
    lowerbound,upperbound = extrema(_log_rand_unif(splits[i-1],splits[i]))
    #subsample 
    push!(inds,lowerbound-1+argmin(data[lowerbound:upperbound]))
    # push!(bounds,(lowerbound,upperbound))
  end
  #do last chunk
  if maxlen>splits[end] #didn't already do this one
    # lowerbound,upperbound = extrema(rand(splits[end]:maxlen,2))
    lowerbound,upperbound = extrema(_log_rand_unif(splits[end],maxlen))
    
    push!(inds,lowerbound-1+argmin(data[lowerbound:upperbound]))
    # push!(bounds,(lowerbound,upperbound))
  end
  return inds#,bounds
end

"""
    epidemic_sample3(E::Union{EventEpidemic.SEIRData,EventEpidemic.SIRData},
                    ktrials::Int,
                    kfailures::Int)

runs ktrials of an epidemic from the same seed node and return the number of 
times each node was infected.
"""
function epidemic_sample3(E::Union{EventEpidemic.SEIRData,EventEpidemic.SIRData},
                    ktrials::Int,
                    kfailures::Int)

  #run ntrials and if we have more than kfailures, pick a different seed node and try again
  seed_node = rand(1:length(E.snodes))
  z = zeros(Float64,length(E.snodes))
  curr_failures = 0
  for trial_num = 1:ktrials
    l,E = EventEpidemic.epidemic(E,seed_node)
    netinfs, newinfs = EventEpidemic.get_infdata(l, E)
    while  sum(newinfs) < min(500,0.01*lastindex(z)) && curr_failures<kfailures
      curr_failures +=1 
      l,E = EventEpidemic.epidemic(E,seed_node)
      netinfs, newinfs = EventEpidemic.get_infdata(l, E)
    end

    if curr_failures<kfailures
        z.+=min.(E.itime,l+1)
    end
  end
  
  
  #figure out sampling for this graph 
  if curr_failures<kfailures 
    #rank nodes, take sweepcut, and subsample
    p = sweepcut(A,-z)
    inds = _subsample2(p.conductance)

    #get stats and push to thread-local storage
    for ind in inds
      push!(ncpdata,[seed_node,ind,p.conductance[ind],p.cut[ind],p.volume[ind]])
    end
  end
end

function tmpfcn3()
  #update epidemic parameters 
  update_E!(E,0.3,0.05,0)
  epidemic_sample3(E,20,5)   
end

function saving_data3(gname::String,dst::String,worker::Int,model::String="seir")
  #saves ncpdata to files 
  @spawnat worker CSV.write(joinpath(dst,"ncpinfo-epidemic-subsampling-4-$model-$(gname[1:end-5])-worker-$(worker).txt"),ncpdata)
  return true
end

function new_parallel_ncp_code3(gname::String,
                                total_trials::Int=2000;model::String="seir",dst::String="")

    #send data to all threads
    @everywhere ncpdata = DataFrame(seed = Int64[], size=Int64[],
                    cond = Float64[], cut = Float64[], volume_seed = Float64[])
    @everywhere gname = $gname
    @everywhere A = loadGraph($gname,"pipeline/graphs/",sym=true)
    if model in ["seir","sir"]
        @everywhere E = EventEpidemic.SEIRData(A,tmax=100000,beta=0.3) #stronger epidemic should be fine
    else
        error("model should be one of sir or seir")
    end
    

    #define parameters
    ps = 1:total_trials
    #run diffusions
    @showprogress pmap(x->tmpfcn3(),ps)

    #save data (from each thread)
    println("diffusions completed. saving data.")
    pmap(x->saving_data3(gname,dst,x,model),workers())
   
    println("data saved, combining files")
    #combine files 
    
    fnames = filter(x->startswith(x,"ncpinfo-epidemic-subsampling-4-$model-$(gname[1:end-5])-worker"),readdir(dst))
    #serially load files and append to main files
    ncp,headerinfo = readdlm(joinpath(dst,"$(fnames[1])"),',',header=true)
    for f in fnames[2:end]
        tmp,headerinfo = readdlm(joinpath(dst,"$(f)"),',',header=true)
        ncp = vcat(ncp,tmp)
    end
    ncp = DataFrame(ncp,vec(headerinfo))
    CSV.write(joinpath(dst,"ncpinfo-epidemic-subsampling-4-$model-$total_trials-$(gname[1:end-5]).txt"),ncp)
    #remove the worker files 
    rm.(joinpath.(dst,fnames))
    return true
end 
