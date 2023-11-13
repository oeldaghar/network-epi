#parallel_epidemic_ncp
#initialize graph and epidemic data structure
#pmap(beta,trial) -> run diffusion and return relevant results
#assumes A and E are loaded locally on processor
using Distributed, ProgressMeter, CSV

using MatrixNetworks, LinearAlgebra, SparseArrays, DataFrames

using DelimitedFiles
using CSV 

# function epidemic_ncp_sample(E::Union{EventEpidemic.SIRData,EventEpidemic.SEIRData};getsets::Bool=false)
#     seed_node = rand(1:size(A,1))

#     #nontrivial epidemic
#     l,E = EventEpidemic.epidemic(E,seed_node)
#     net,new = EventEpidemic.get_infdata(l,E)
#     while sum(new)==0
#         l,E = EventEpidemic.epidemic(E,seed_node)
#         net,new = EventEpidemic.get_infdata(l,E)
#     end
#     #do sweepcut
#     p = MatrixNetworks.sweepcut(A,-1*E.itime)
#     ind = argmin(p.conductance)
#     if ind > size(A,1)/2
#         bs = p.p[ind+1:end]
#     else
#         bs = p.p[1:ind]
#     end
#     #return relevant info
#     p = MatrixNetworks.sweepcut(A,-1*E.itime)
#     ind = argmin(p.conductance)
#     slen = Int(length(bs))
#     if getsets==false
#         bs = Vector{Int}()
#     end

#     return bs,(seed_node,E.beta,slen,p.conductance[ind],Float64(p.cut[ind]),p.volume[ind],nnz(A)-p.volume[ind])
# end

# function compute_time(beta::Float64;getsets::Bool=false)
#     # t1 = time() #@elapsed messes with the scoping of variables
#     dt = @elapsed update_E!(E,beta,E.gamma,0)
#     dt += @elapsed ncpinfo = epidemic_ncp_sample1(E)#,getsets=getsets)
#     # dt = time()-t1
#     return [collect(ncpinfo);dt]#, Set(bs)
# end

# function parallel_epidemic_ncp(A::SparseMatrixCSC,ntrials::Int;getsets::Bool=false)
#     """
#         fucntion for using epidemic results to rank node importance prior to performing a sweepcut
#     """
#     Avol=sum(A)
#     E = EventEpidemic.SEIRData(A)
#     #initialize info
#     ncpdata = DataFrame(seed = Int64[], eps = Float64[], size=Int64[],
#                     cond = Float64[], cut = Float64[], volume_seed = Float64[],
#                     volume_other = Float64[], ComputeTime = Float64[])
#     bestsets = Vector{Set{Int64}}()
#     @everywhere @eval A = $A
#     @everywhere @eval E = EventEpidemic.SEIRData(A)
#     @everywhere @eval Avol=$Avol

#     betas = [1e-2 3e-2 5e-2 8e-2 1e-1 5e-1]
#     ps = vcat(Base.Iterators.product(betas,1:ntrials)...)

#     # @everywhere tfcn1((beta,ind)) = compute_time(beta)
#     #record info
#     data = @showprogress pmap(c->compute_time(c[1],getsets=getsets),ps)
#     push!(ncpdata,first.(data)...)
#     push!(bestsets,map(x->x[2],data)...)
#     return ncpdata,bestsets
# end

# function make_buckets(arr::Vector{Int},bounds::Vector{Int})
#     """split arr according to bounds. so i in bucket k if arr[i]>=bounds[k] but arr[i]<bounds[k+1]
#     this sorts bounds first.
#     remaining nodes placed in last bucket. #TODO: handle edge case where arr[i]<bounds[1]
#     """
#     #split arr into buckets 
#     sort!(bounds)
#     bins = Dict{Int,Vector{Int}}()
#     for i = 1:length(bounds)
#         bins[i] = Vector{Int}()
#     end
#     nbins = length(bounds)

#     for (i,v) in enumerate(arr)
#         t=1
#         while t<nbins && v>bounds[t]
#           t+=1
#         end
#         push!(bins[t],i)
#     end
#     return bins
# end
# make_buckets(arr::Vector{Int},bounds::Vector{Float64}) = make_buckets(arr,Int.(bounds))

# function epidemic_ncp_sample1(E::Union{EventEpidemic.SIRData,EventEpidemic.SEIRData})
#     """perform a single epidemic diffusion and samples sets using sweepcut profile of MatrixNetworks"""
#     seed_node = rand(1:size(A,1))

#     #nontrivial epidemic
#     z = zeros(Int,length(E.snodes))
#     for i=1:10
#         l,E = EventEpidemic.epidemic(E,seed_node,tmax=25000)
#         net,new = EventEpidemic.get_infdata(l,E)
#         while sum(new)==0 #conditioning the epidemic on infecting at least one other person
#             l,E = EventEpidemic.epidemic(E,seed_node)
#             net,new = EventEpidemic.get_infdata(l,E)
#         end

#         z[E.ctime .<= l] .+= 1
#     end
#     #---------- modify code here to change how we sample sets from a diffusion -----------#
#     #get conductance for sets from sweepcut
    
#     # sx = -1*E.itime #MatrixNetworks uses node ranking 
#     # p = MatrixNetworks.sweepcut(A,sx)
#     p = MatrixNetworks.sweepcut(A,z)
   
#     #using chain  
#     # px = sortperm(p.conductance)
#     #splitting nodes into buckets for sampling 
#     # pbins = make_buckets(px[1:end-1],[1;10;1e2;1e3;1e4;1e5])
#     # sizes = map(x->px[x],values(pbins)) 
#     # pss = map(x->p.conductance[sizes[x]],1:length(sizes))
#     #get best conductance set in each bucket 
#     # inds = map(x->sizes[x][argmin(pss[x])],1:length(sizes))

#     inds = argmin(p)
#     #---------- modify code here to change what we record. don't recommend changing this. -----------#
#     #make relevant ncpinfo
#     result = Vector{Tuple}()
#     for ind in inds
#         push!(result, (seed_node,E.beta,ind,p.conductance[ind],Float64(p.cut[ind]),p.volume[ind],nnz(A)-p.volume[ind]))
#     end
#     return result 
# end

# function compute_time1(beta_timedelay::Tuple{Float64,Float64})
#     """
#         updates epidemic parameters for E, records time to do a single diffusion, and returns ncpinfo 
#     """
#     beta,a = beta_timedelay[1],beta_timedelay[2]
#     dt = @elapsed update_E!(E,beta,E.gamma,0)
#     #update_E!(E,beta,gamma,quarantine_capacity) #dont change quarantine_capacity. 
#     dt += @elapsed E.a = a
#     dt += @elapsed ncpinfo = epidemic_ncp_sample1(E)
#     return push!.(collect.(ncpinfo),dt)
# end

# function parallel_epidemic_ncp1(A::SparseMatrixCSC,ntrials::Int)
#     """
#         fucntion for using epidemic results to rank node importance prior to performing a sweepcut
#     """
#     Avol=sum(A)
#     E = EventEpidemic.SEIRData(A,tmax=100000)
#     #initialize info
#     ncpdata = DataFrame(seed = Int64[], eps = Float64[], size=Int64[],
#                     cond = Float64[], cut = Float64[], volume_seed = Float64[],
#                     volume_other = Float64[], ComputeTime = Float64[])
#     #send data to all threads
#     @everywhere @eval A = $A
#     @everywhere @eval E = EventEpidemic.SEIRData(A,tmax=100000)
#     @everywhere @eval Avol=$Avol

#     #---------------change code here to change parameters for epidemics-------------------
#     #changing beta 
#     betas = [5e-3 1e-2 3e-2] #infection probabilities
#     as = [0.1 10.0 500.0] #parameter for exposed->infected time delay. larger values indicate larger delays in expectation
    
#     ps = vcat(Base.Iterators.product(betas,as)...)
#     ps = vcat(Base.Iterators.product(ps,1:ntrials)...)
#     #ps = [ ( (beta,a) ,trial_number) ]

#     #main function for parallelism 
#     data = @showprogress pmap(c->compute_time1(c[1]),ps)
    
#     #flattening data
#     push!(ncpdata,vcat(data...)...)
#     return ncpdata
# end


## plotting fcns
# include("hexbins.jl")
# include("ncpplots.jl")

# function plot_hexbinncp(ncp::DataFrame,A::SparseMatrixCSC)
#     f = myncpplot(map(x->min(ncp.size[x],size(A,1)-ncp.size[x]),1:size(ncp,1)),ncp.cond)
#     Plots.plot!(xlabel="size",ylabel="conductance",xlims=(5,1.2*size(A,1)),
#         ylims=(1e-5,1),label="")
#
#     f
# end
#
# function make_hexbinncp_plots1(gname::String,ntrials::Int;gpath::String="")
#     """
#         perform ncp diffusions and generate hexbin plot  
#     """
#     A = loadGraph(gname,gpath)
#     davg = nnz(A)/size(A,1)
#     println("doing ncp calculation")
#     ncp = parallel_epidemic_ncp1(A,ntrials);
    
#     println("making hexbin ncp plot")
#     f = plot_hexbinncp(ncp,A);
#     title!(f,"NCP - $(gname[1:end-5])\ndavg: $(round(davg,digits=2)), nnodes: $(size(A,1))")
#     # Plots.savefig(joinpath(dst,"ncphexbin-epidemic1-$ntrials-$(gname[1:end-5]).png"))

#     #save raw data
#     # println("saving ncpinfo")
#     # CSV.write(joinpath(dst,"ncpinfo-epidemic1-$ntrials-$(gname[1:end-5]).txt"),ncp)
#     return ncp,f
# end




#epidemic
#sweepcut
#append data
#save to file (using threadid to avoid fname issues)
function epidemic_sample(E::Union{EventEpidemic.SEIRData,EventEpidemic.SIRData},
                    target_infs::Int,
                    ktrials::Int,
                    kfailures::Int)
  """
    runs ktrials of an epidemic from the same seed node and return the number of 
    times each node was infected.
  """

  #run ntrials and if we have more than kfailures, pick a different seed node and try again
  seed_node = rand(1:length(E.snodes))
  z = zeros(Int,length(E.snodes))
  curr_failures = 0
  for trial_num = 1:ktrials
    # println("working on trials $trial_num")
    l,E = EventEpidemic.partial_epidemic(E,seed_node,target_infs)
    netinfs, newinfs = EventEpidemic.get_infdata(l, E)
    while  sum(newinfs) < target_infs && curr_failures<kfailures
      curr_failures +=1 
      l,E = EventEpidemic.partial_epidemic(E,seed_node,target_infs)
      netinfs, newinfs = EventEpidemic.get_infdata(l, E)
    end
    if curr_failures<kfailures
        z.+=E.snodes  
    end
  end
  #push data to thread-local storage
  #figure out sampling for this graph 
  maxlen = size(A,1)-1
  sampling_splits = [1;100;1000;10000]
  for i=1:lastindex(sampling_splits)
    if sampling_splits[i]>maxlen
      sampling_splits = sampling_splits[1:i-1]
      break
    end
  end
  if ktrials-curr_failures>=10 
    p = sweepcut(A, ktrials .- z);
    #forced subsampling based on size of set. 
    inds = Vector{Int}()
    for i=2:lastindex(sampling_splits)
      push!(inds,sampling_splits[i-1]-1+argmin(p.conductance[sampling_splits[i-1]:sampling_splits[i]]))
    end
    #do tail end of graph
    if maxlen>sampling_splits[end] #didn't already do this one
      push!(inds,sampling_splits[end]-1+argmin(p.conductance[sampling_splits[end]:end]))
    end
    # if maxlen>=100
    #   push!(inds,argmin(p.conductance[1:100]))
    # else
    #   push!(inds,argmin(p.conductance[1:maxlen]))
    # end
    # if maxlen>=1000
    # push!(inds,99+argmin(p.conductance[100:1000]))
    # #quick fix. need to do this in a better way
    # if length(p.conductance)>=10000
    #   push!(inds,999+argmin(p.conductance[1000:10000]))
    #   push!(inds,9999+argmin(p.conductance[10000:end]))
    # else
    #   push!(inds,999+argmin(p.conductance[1000:end]))
    # end
    for ind in inds
        push!(ncpdata,[seed_node,ind,p.conductance[ind],p.cut[ind],p.volume[ind], target_infs])
    end
  end
  return ktrials-curr_failures>=10  
end

function epidemic_sample1(E::Union{EventEpidemic.SEIRData,EventEpidemic.SIRData},
                    target_infs::Int,
                    ktrials::Int,
                    kfailures::Int)
  """
    runs ktrials of an epidemic from the same seed node and return the number of 
    times each node was infected.
  """

  #run ntrials and if we have more than kfailures, pick a different seed node and try again
  seed_node = rand(1:length(E.snodes))
  z = zeros(Int,length(E.snodes))
  curr_failures = 0
  for trial_num = 1:ktrials
    # println("working on trials $trial_num")
    l,E = EventEpidemic.partial_epidemic(E,seed_node,target_infs)
    netinfs, newinfs = EventEpidemic.get_infdata(l, E)
    while  sum(newinfs) < target_infs && curr_failures<kfailures
      curr_failures +=1 
      l,E = EventEpidemic.partial_epidemic(E,seed_node,target_infs)
      netinfs, newinfs = EventEpidemic.get_infdata(l, E)
    end
    if curr_failures<kfailures
        # z.+=E.snodes  
        z.+=min.(E.itime,l+1)
    end
  end
  #smooth out the tail end of z 

  #push data to thread-local storage
  #figure out sampling for this graph 
  maxlen = size(A,1)-1
  sampling_splits = [1;100;1000;10000]
  for i=1:lastindex(sampling_splits)
    if sampling_splits[i]>maxlen
      sampling_splits = sampling_splits[1:i-1]
      break
    end
  end
  if ktrials-curr_failures>=10 
    p = sweepcut(A, ktrials .- z);
    #forced subsampling based on size of set. 
    #can get accumulation of sets at the boundary of these cutoffs (think spatial network)
    #to counter this, do the subsampling by 
    inds = Vector{Int}()
    for i=2:lastindex(sampling_splits)
      #subsample upperbound
      upperbound = rand(sampling_splits[i-1]:sampling_splits[i])
      #subsample 
      push!(inds,sampling_splits[i-1]-1+argmin(p.conductance[sampling_splits[i-1]:upperbound]))
    end
    #do tail end of graph
    if maxlen>sampling_splits[end] #didn't already do this one
      upperbound = rand(sampling_splits[end]:maxlen)
      push!(inds,sampling_splits[end]-1+argmin(p.conductance[sampling_splits[end]:upperbound]))
    end
    #get stats 
    for ind in inds
        push!(ncpdata,[seed_node,ind,p.conductance[ind],p.cut[ind],p.volume[ind], target_infs])
    end
  end
  return ktrials-curr_failures>=10  
end

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

function epidemic_sample2(E::Union{EventEpidemic.SEIRData,EventEpidemic.SIRData},
                    target_infs::Int,
                    ktrials::Int,
                    kfailures::Int)
  """
    runs ktrials of an epidemic from the same seed node and return the number of 
    times each node was infected.
  """

  #run ntrials and if we have more than kfailures, pick a different seed node and try again
  seed_node = rand(1:length(E.snodes))
  z = zeros(Int,length(E.snodes))
  curr_failures = 0
  for trial_num = 1:ktrials
    # println("working on trials $trial_num")
    l,E = EventEpidemic.partial_epidemic(E,seed_node,target_infs)
    netinfs, newinfs = EventEpidemic.get_infdata(l, E)
    while  sum(newinfs) < target_infs && curr_failures<kfailures
      curr_failures +=1 
      l,E = EventEpidemic.partial_epidemic(E,seed_node,target_infs)
      netinfs, newinfs = EventEpidemic.get_infdata(l, E)
    end
    if curr_failures<kfailures
        # z.+=E.snodes  
        z.+=min.(E.itime,l+1)
    end
  end
  
  #push data to thread-local storage
  #figure out sampling for this graph 
  if ktrials-curr_failures>=10 
    p = sweepcut(A, ktrials .- z);

    inds = _subsample2(p.conductance)
    #get stats 
    for ind in inds
        push!(ncpdata,[seed_node,ind,p.conductance[ind],p.cut[ind],p.volume[ind], target_infs])
    end
  end
  return ktrials-curr_failures>=10  
end


function tmpfcn(target_infs::Int)
  #update epidemic parameters 
  update_E!(E,0.2,0.05,0)
  epidemic_sample(E,target_infs,30,5)   
end

function tmpfcn1(target_infs::Int)
  #update epidemic parameters 
  update_E!(E,0.2,0.05,0)
  epidemic_sample1(E,target_infs,20,5)   
end

function tmpfcn2(target_infs::Int)
  #update epidemic parameters 
  update_E!(E,0.2,0.05,0)
  epidemic_sample2(E,target_infs,20,5)   
end


function saving_data(gname::String,dst::String,worker::Int,model::String="seir")
  #saves ncpdata to files 
  @spawnat worker CSV.write(joinpath(dst,"ncpinfo-epidemic-$model-$(gname[1:end-5])-worker-$(worker).txt"),ncpdata)
  return true
end

function saving_data1(gname::String,dst::String,worker::Int,model::String="seir")
  #saves ncpdata to files 
  @spawnat worker CSV.write(joinpath(dst,"ncpinfo-epidemic-subsampling-$model-$(gname[1:end-5])-worker-$(worker).txt"),ncpdata)
  return true
end

function saving_data2(gname::String,dst::String,worker::Int,model::String="seir")
  #saves ncpdata to files 
  @spawnat worker CSV.write(joinpath(dst,"ncpinfo-epidemic-subsampling-2-$model-$(gname[1:end-5])-worker-$(worker).txt"),ncpdata)
  return true
end

function new_parallel_ncp_code(gname::String,
                                total_trials::Int=2000;model::String="seir",dst::String="")

    
    #send data to all threads
    @everywhere ncpdata = DataFrame(seed = Int64[], size=Int64[],
                    cond = Float64[], cut = Float64[], volume_seed = Float64[], target = Int64[])
    @everywhere gname = $gname
    @everywhere A = loadGraph($gname,"pipeline/graphs/",sym=true)
    if model=="seir"
        @everywhere E = EventEpidemic.SEIRData(A,tmax=100000,beta=0.1) #stronger epidemic should be fine
    elseif model == "sir"
        @everywhere E = EventEpidemic.SIRData(A,tmax=100000,beta=0.5) #stronger epidemic should be fine
    else
        error("model should be one of sir or seir")
    end
    

    #define parameters
    # target_infs = [10;50;100;200;500;1000;2000;5000;10000]
    target_infs = [200;500]
    max_ind = floor(Int,total_trials/length(target_infs))
    ps = Base.product(target_infs,1:max_ind)
    #run diffusions
    #@showprogress pmap(x->...,parameters) 
    @showprogress pmap(x->tmpfcn(x[1]),ps)

    #save data (from each thread)
    println("diffusions completed. saving data.")
    pmap(x->saving_data(gname,dst,x,model),workers())
   
    println("data saved, combining files")
    #combine files 
    
    fnames = filter(x->startswith(x,"ncpinfo-epidemic-$model-$(gname[1:end-5])-worker"),readdir(dst))
    #serially load files and append to main files
    # using DelimitedFiles
    ncp,headerinfo = readdlm(joinpath(dst,"$(fnames[1])"),',',header=true)
    for f in fnames[2:end]
        tmp,headerinfo = readdlm(joinpath(dst,"$(f)"),',',header=true)
        ncp = vcat(ncp,tmp)
    end
    ncp = DataFrame(ncp,vec(headerinfo))
    CSV.write(joinpath(dst,"ncpinfo-epidemic-$model-$total_trials-$(gname[1:end-5]).txt"),ncp)
    #remove the worker files 
    rm.(joinpath.(dst,fnames))
end 


function new_parallel_ncp_code1(gname::String,
                                total_trials::Int=2000;model::String="seir",dst::String="")

    
    #send data to all threads
    @everywhere ncpdata = DataFrame(seed = Int64[], size=Int64[],
                    cond = Float64[], cut = Float64[], volume_seed = Float64[], target = Int64[])
    @everywhere gname = $gname
    @everywhere A = loadGraph($gname,"pipeline/graphs/",sym=true)
    if model=="seir"
        @everywhere E = EventEpidemic.SEIRData(A,tmax=100000,beta=0.5) #stronger epidemic should be fine
    elseif model == "sir"
        @everywhere E = EventEpidemic.SIRData(A,tmax=100000,beta=0.5) #stronger epidemic should be fine
    else
        error("model should be one of sir or seir")
    end
    

    #define parameters
    target_infs = [200;500]
    max_ind = floor(Int,total_trials/length(target_infs))
    ps = Base.product(target_infs,1:max_ind)
    #run diffusions
    #@showprogress pmap(x->...,parameters) 
    @showprogress pmap(x->tmpfcn1(x[1]),ps)

    #save data (from each thread)
    println("diffusions completed. saving data.")
    pmap(x->saving_data1(gname,dst,x,model),workers())
   
    println("data saved, combining files")
    #combine files 

    fnames = filter(x->startswith(x,"ncpinfo-epidemic-subsampling-$model-$(gname[1:end-5])-worker"),readdir(dst))
    #serially load files and append to main 
    # using DelimitedFiles
    ncp,headerinfo = readdlm(joinpath(dst,"$(fnames[1])"),',',header=true)
    for f in fnames[2:end]
        tmp,headerinfo = readdlm(joinpath(dst,"$(f)"),',',header=true)
        ncp = vcat(ncp,tmp)
    end
    ncp = DataFrame(ncp,vec(headerinfo))
    CSV.write(joinpath(dst,"ncpinfo-epidemic-subsampling-$model-$total_trials-$(gname[1:end-5]).txt"),ncp)
    #remove the worker files 
    rm.(joinpath.(dst,fnames))
end 


function new_parallel_ncp_code2(gname::String,
                                total_trials::Int=2000;model::String="seir",dst::String="")

    
    #send data to all threads
    @everywhere ncpdata = DataFrame(seed = Int64[], size=Int64[],
                    cond = Float64[], cut = Float64[], volume_seed = Float64[], target = Int64[])
    @everywhere gname = $gname
    @everywhere A = loadGraph($gname,"pipeline/graphs/",sym=true)
    if model=="seir"
        @everywhere E = EventEpidemic.SEIRData(A,tmax=100000,beta=0.5) #stronger epidemic should be fine
    elseif model == "sir"
        @everywhere E = EventEpidemic.SIRData(A,tmax=100000,beta=0.5) #stronger epidemic should be fine
    else
        error("model should be one of sir or seir")
    end
    

    #define parameters
    target_infs = [200;500]
    max_ind = floor(Int,total_trials/length(target_infs))
    ps = Base.product(target_infs,1:max_ind)
    #run diffusions
    #@showprogress pmap(x->...,parameters) 
    @showprogress pmap(x->tmpfcn2(x[1]),ps)

    #save data (from each thread)
    println("diffusions completed. saving data.")
    pmap(x->saving_data2(gname,dst,x,model),workers())
   
    println("data saved, combining files")
    #combine files 

    fnames = filter(x->startswith(x,"ncpinfo-epidemic-subsampling-2-$model-$(gname[1:end-5])-worker"),readdir(dst))
    #serially load files and append to main 
    # using DelimitedFiles
    ncp,headerinfo = readdlm(joinpath(dst,"$(fnames[1])"),',',header=true)
    for f in fnames[2:end]
        tmp,headerinfo = readdlm(joinpath(dst,"$(f)"),',',header=true)
        ncp = vcat(ncp,tmp)
    end
    ncp = DataFrame(ncp,vec(headerinfo))
    CSV.write(joinpath(dst,"ncpinfo-epidemic-subsampling-2-$model-$total_trials-$(gname[1:end-5]).txt"),ncp)
    #remove the worker files 
    rm.(joinpath.(dst,fnames))
end 

### 3rd attempt
function _subsample3(data,splits=[1;100;1000;10000;100000])
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
    tmpbins = collect(splits[i-1]:(splits[i]-splits[i-1])/100:splits[i])
    tmpbins = unique(floor.(Int,tmpbins))
    ind = rand(1:99)
    lowerbound,upperbound = tmpbins[ind],tmpbins[ind+1]
    #subsample 
    push!(inds,lowerbound-1+argmin(data[lowerbound:upperbound]))
    # push!(bounds,(lowerbound,upperbound))
  end
  #do last chunk
  if maxlen>splits[end] #didn't already do this one

    tmpbins = collect(splits[end]:(maxlen-splits[end])/100:maxlen)
    tmpbins = unique(floor.(Int,tmpbins))
    ind = rand(2:min(100,length(tmpbins)-1))

    lowerbound,upperbound = tmpbins[ind],maxlen
    #subsample     
    push!(inds,lowerbound-1+argmin(data[lowerbound:upperbound]))
    # push!(bounds,(lowerbound,upperbound))
  end
  return inds#,bounds
end

function epidemic_sample3(E::Union{EventEpidemic.SEIRData,EventEpidemic.SIRData},
                    ktrials::Int,
                    kfailures::Int)
  """
    runs ktrials of an epidemic from the same seed node and return the number of 
    times each node was infected.
  """

  #run ntrials and if we have more than kfailures, pick a different seed node and try again
  seed_node = rand(1:length(E.snodes))
  z = zeros(Float64,length(E.snodes))
  curr_failures = 0
  for trial_num = 1:ktrials
    # println("working on trials $trial_num")
    # l,E = EventEpidemic.partial_epidemic(E,seed_node,target_infs)
    l,E = EventEpidemic.epidemic(E,seed_node)
    netinfs, newinfs = EventEpidemic.get_infdata(l, E)
    while  sum(newinfs) < min(500,0.01*lastindex(z)) && curr_failures<kfailures
      curr_failures +=1 
      # l,E = EventEpidemic.partial_epidemic(E,seed_node,target_infs)
      l,E = EventEpidemic.epidemic(E,seed_node)
      netinfs, newinfs = EventEpidemic.get_infdata(l, E)
    end
    #define weight fcn for nodes 
    weight_fcn(x) = E.snodes[x] ? 0.0 : 1/(1+E.itime[x])
    if curr_failures<kfailures
        # z.+=E.snodes  
        z.+=min.(E.itime,l+1)
        # z .+= weight_fcn.(1:lastindex(E.snodes))
    end
  end
  
  #push data to thread-local storage
  #figure out sampling for this graph 
  # if ktrials-curr_failures>=10 
  if curr_failures<kfailures 
    #smooth out epidemic to rest of graph using mean field approach
    # w = (ktrials .- z)./ktrials
    # p = sweepcut(A, ktrials .- z);

    # p = sweepcut(A,z)
    p = sweepcut(A,-z)
    inds = _subsample2(p.conductance)
    # inds = _subsample3(p.conductance)

    #get stats 
    for ind in inds
      push!(ncpdata,[seed_node,ind,p.conductance[ind],p.cut[ind],p.volume[ind]])
    end
  end
  # return ktrials-curr_failures>=  
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
    # target_infs = [10;50;100;200;500;1000;2000;5000;10000]
    # target_infs = [200;500]
    # max_ind = floor(Int,total_trials/length(target_infs))
    ps = 1:total_trials
    #run diffusions
    #@showprogress pmap(x->...,parameters) 
    @showprogress pmap(x->tmpfcn3(),ps)

    #save data (from each thread)
    println("diffusions completed. saving data.")
    pmap(x->saving_data3(gname,dst,x,model),workers())
   
    println("data saved, combining files")
    #combine files 
    
    fnames = filter(x->startswith(x,"ncpinfo-epidemic-subsampling-4-$model-$(gname[1:end-5])-worker"),readdir(dst))
    #serially load files and append to main files
    # using DelimitedFiles
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
