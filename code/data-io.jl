using ProgressMeter
# using SparseArrays

## functions for loading in infection data 
#for loading in data 
function read_inf_data(gname;dloc=dloc,beta=0.1,gamma=0.05,method="sir",dtype="cinfs",exp=5.0)
    if endswith(gname,".smat")
        g = gname[1:end-5]
    else
        g = gname
    end

    if beta!=0
        fname = dloc*"$dtype-$g-$method-$beta-$gamma.txt"
    else
        fname = dloc*"$dtype-$g-$method-$exp.txt"
    end
    
    if lowercase(dtype)=="cinfs" || lowercase(dtype)=="tinfs"
        res = Vector{Vector{Int}}()
        open(fname,"r") do io
            while !eof(io)
                push!(res,parse.(Int,split(readline(io),",")))
            end
        end
    else
        error("no such file or wrong data format. this reads Vector{Vector{Int}} data")
    end
    return res
end

## functions for getting parameters from diffusion data 

function sort_fnames(gname::String;
                beta::Float64=0.05,
                gamma::Float64=0.05,
                # rps::Vector{Float64}=vec([0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.25 0.5 1.0 10.0 100.0]),
                dloc::String="pipeline/data/",
                dtype::String="tinfs",
                method::String="seir",
                rewiring_type::String="-rewired-",
                diffusion_type::String="uniform")
    """
        given gname and diffusion parameters, get graphs of that type and sort the file names by rewiring fraction

        returns rewired_gnames and rewiring_fractions
    """


    if endswith(gname,".smat")
        gname = gname[1:end-5]
    end
    rewiring_type = strip(rewiring_type,'-')

    dataDir = joinpath(dloc,"$gname/diffusions/$diffusion_type/")

    #sort gnames by rewiring parameter
    gnames = readdir(dataDir)
    filter!(x->startswith(x,"$dtype-$rewiring_type-"),gnames) #excludes original graph     
    filter!(x->occursin("$method-$beta-$gamma",x),gnames)
    ind = length("$dtype-$rewiring_type-")+1

    rewiring_fractions = map(x->first(split(x[ind:end],"-")),gnames)
    p = sortperm(parse.(Float64,rewiring_fractions))

    return gnames[p],Vector{String}(rewiring_fractions[p])
end
#get all betas associated with a graph type 
#currently hard-coded for 
# - "uniform" diffusions 
# - "seir" model

function get_betas(gname::String)
    g = canonical_graph_name(gname)[1:end-5]
    #find all betas
    tmp_fnames = filter(x->startswith(x,"tinfs-$g-seir-"),readdir("pipeline/data/$g/diffusions/uniform/"))
    tmp_betas = map(x->split(x,"tinfs-$g-seir-")[2],tmp_fnames)

    #extract betas
    tmp_betas = map(x->split(x,"-")[1],tmp_betas)
    betas = parse.(Float64,tmp_betas)
    return betas
end

##

function load_single_rewiring_data(gname::String;
                            beta::Float64=0.05,
                            gamma::Float64=0.05,
                            # rps::Vector{Float64}=vec([0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.25 0.5 1.0 10.0 100.0]),
                            dloc::String="pipeline/data/",
                            dtype::String="tinfs",
                            method::String="seir",
                            rewiring_type::String="triange-rewired",
                            diffusion_type::String="triangle-weighted")
    """
        given gname and other diffusion paramaters, load in an array X of the format

        rewiring_percent (rows) vs (qpercent \times nnodes in trial) (cols)
    """

    if endswith(gname,".smat")
        gname = gname[1:end-5]
    end

    dataDir = joinpath(dloc,"$gname/diffusions/$diffusion_type/")
    rewiring_type = String(strip(rewiring_type,'-'))
    
    gnames,rewiring_fractions = sort_fnames(gname,beta=beta,gamma=gamma,dloc=dloc,dtype=dtype,method=method,
                                        rewiring_type=rewiring_type,diffusion_type=diffusion_type)

    data = Vector{Vector{Vector{Int}}}()

    for i in eachindex(gnames)
        push!(data,cumsum.(read_inf_data("$rewiring_type-$(rewiring_fractions[i])-$gname",dloc=dataDir,beta=beta,gamma=gamma,dtype=dtype,method=method)))
    end
    return data,rewiring_fractions
end


function load_double_rewiring_data(gname::String;
                            beta::Float64=0.05,
                            gamma::Float64=0.05,
                            # rps::Vector{Float64}=vec([0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.25 0.5 1.0 10.0 100.0]),
                            dloc::String="pipeline/data/",
                            dtype::String="tinfs",
                            method::String="seir",
                            rewiring_type1::String="rewired",
                            rewiring_type2::String="er",
                            diffusion_type::String="uniform")
    """
        given gname and other diffusion paramaters, load in an array X of the format

        rewiring_percent (rows) vs (qpercent \times nnodes in trial) (cols)
    """

    if endswith(gname,".smat")
        gname = gname[1:end-5]
    end
    dataDir = joinpath(dloc,"$gname/diffusions/$diffusion_type/")

    rewiring_type1 = String(strip(rewiring_type1,'-'))
    rewiring_type2= String(strip(rewiring_type2,'-'))

    gnames1,rewiring_fractions1 = sort_fnames(gname,rewiring_type=rewiring_type1,
                                        beta=beta,gamma=gamma,dloc=dloc,dtype=dtype,method=method,diffusion_type=diffusion_type)
    gnames1 = reverse(gnames1)
    rewiring_fractions1 = reverse(rewiring_fractions1)

    gnames2,rewiring_fractions2 = sort_fnames(gname,rewiring_type=rewiring_type2,
                                        beta=beta,gamma=gamma,dloc=dloc,dtype=dtype,method=method,diffusion_type=diffusion_type)

    # gnames = [collect(gnames1);gname;collect(gnames2)]

    data = Vector{Vector{Vector{Int}}}()
    #first rewiring type 
    for i in eachindex(gnames1)
        push!(data,cumsum.(read_inf_data("$rewiring_type1-$(rewiring_fractions1[i])-$gname",dloc=dataDir,beta=beta,gamma=gamma,dtype=dtype,method=method)))
    end
    #do original graph
    push!(data,cumsum.(read_inf_data("$gname",dloc=dataDir,beta=beta,gamma=gamma,dtype=dtype,method=method)))
    #second rewiring type 
    for i in eachindex(gnames2)
        push!(data,cumsum.(read_inf_data("$rewiring_type2-$(rewiring_fractions2[i])-$gname",dloc=dataDir,beta=beta,gamma=gamma,dtype=dtype,method=method)))
    end
    
    return data,rewiring_fractions1,rewiring_fractions2
end


#functions for loading in full and aggregated diffusion data 
#will need to eventually change these function names 
function get_plotting_data(gname::String;
                    beta::Float64=0.05,
                    gamma::Float64=0.05,
                    gpath::String="pipeline/graphs/",
                    dloc::String="pipeline/data/",
                    dtype::String="tinfs",
                    method::String="seir",
                    rewiring_type1::String="rewired",
                    rewiring_type2::String="er",
                    diffusion_type::String="uniform",
                    ntrials::Int=50,
                    qpercents::Vector{Int}=collect(0:15))
    
    A = loadGraph(gname,gpath)
    nnodes = size(A,1) 
    data,rewiring_fractions1,rewiring_fractions2 = load_double_rewiring_data(gname,
                    beta=beta,gamma=gamma,dloc=dloc,dtype=dtype,method=method,
                    rewiring_type1=rewiring_type1,
                    rewiring_type2=rewiring_type2,
                    diffusion_type=diffusion_type)

    if dtype == "tinfs"
        clabel = "Fraction of Infected Nodes"
        data = hcat(map(x->last.(x),data)...)./nnodes
    elseif dtype == "cinfs"
        clabel = "Normalized Maximum Infection Size"
        data = hcat(map(x->maximum.(x),data)...)./nnodes
    end
    #aggregating over nodes in a cell (same rewiring percent + same quarantine percentage). 
    nrows,ncols = size(data)
    nrows = Int(nrows/ntrials)
    aggreated_data = zeros(Float64,nrows,ncols)    
    for col in eachindex(1:ncols)
        for row in eachindex(1:nrows)
            offset = (row-1)*ntrials
            aggreated_data[row,col] = sum(data[offset+1:offset+ntrials,col])/ntrials
        end
    end
    return data,aggreated_data,clabel,rewiring_fractions1,rewiring_fractions2
end

function get_single_plotting_data(gname::String;
                    beta::Float64=0.05,
                    gamma::Float64=0.05,
                    gpath::String="pipeline/graphs/",
                    dloc::String="pipeline/data/",
                    dtype::String="tinfs",
                    method::String="seir",
                    rewiring_type::String="triangle-rewired",
                    diffusion_type::String="triangle-weighted",
                    ntrials::Int=50,
                    qpercents::Vector{Int}=collect(0:15))
    
    A = loadGraph(gname,gpath)
    nnodes = size(A,1) 
    data,rewiring_fractions = load_single_rewiring_data(gname,
                    beta=beta,gamma=gamma,dloc=dloc,dtype=dtype,method=method,
                    rewiring_type=rewiring_type,
                    diffusion_type=diffusion_type)

    if dtype == "tinfs"
        clabel = "Fraction of Infected Nodes"
        data = hcat(map(x->last.(x),data)...)./nnodes
    elseif dtype == "cinfs"
        clabel = "Normalized Maximum Infection Size"
        data = hcat(map(x->maximum.(x),data)...)./nnodes
    end
    #aggregating over nodes in a cell (same rewiring percent + same quarantine percentage). 
    
    nrows,ncols = size(data)
    nrows = Int(nrows/ntrials)
    aggreated_data = zeros(Float64,nrows,ncols)    
    for col in eachindex(1:ncols)
        for row in eachindex(1:nrows)
            offset = (row-1)*ntrials
            aggreated_data[row,col] = sum(data[offset+1:offset+ntrials,col])/ntrials
        end
    end
    return data,aggreated_data,clabel,rewiring_fractions#,rewiring_fractions2
end


"""
    aggregate_diffusion_data(data,ntrials::Int=50)

aggregates diffusion data to return a matrix of size (# quarantine parameters × # network rewirings)
"""
function aggregate_diffusion_data(data,ntrials::Int=50)
    nrows,ncols = size(data)
    nrows = Int(nrows/ntrials)
    aggreated_data = zeros(Float64,nrows,ncols)    
    for col in eachindex(1:ncols)
        for row in eachindex(1:nrows)
            offset = (row-1)*ntrials
            aggreated_data[row,col] = sum(data[offset+1:offset+ntrials,col])/ntrials
        end
    end
    return aggreated_data
end



"""
    load_triangle_diffusion_data(gname::String
                            beta::Float64=0.05,
                            gamma::Float64=0.05,
                            dloc::String="pipeline/data/",
                            dtype::String="tinfs",
                            method::String="seir",
                            rewiring_type::String="triange-rewired",
                            diffusion_type::String="triangle-weighted")

TBW
"""
function load_triangle_diffusion_data(gname::String,
                            beta::Float64=0.05,
                            gamma::Float64=0.05,
                            dloc::String="pipeline/data/",
                            dtype::String="tinfs",
                            method::String="seir")

    if endswith(gname,".smat")
        gname = gname[1:end-5]
    end

    rewiring_type = "triangle-rewired"
    diffusion_type = "new-triangle-weighted"

    dataDir = joinpath(dloc,"$gname/diffusions/$diffusion_type/")
        
    gnames,rewiring_fractions = sort_fnames(gname,beta=beta,gamma=gamma,dloc=dloc,dtype=dtype,method=method,
                                        rewiring_type=rewiring_type,diffusion_type=diffusion_type)

    data = Vector{Vector{Vector{Int}}}()

    #do original graph
    push!(data,cumsum.(read_inf_data("$gname",dloc=dataDir,beta=beta,gamma=gamma,dtype=dtype,method=method)))
    #do triangle rewired graphs 
    for i in eachindex(gnames)
        push!(data,cumsum.(read_inf_data("$rewiring_type-$(rewiring_fractions[i])-7-$gname",dloc=dataDir,beta=beta,gamma=gamma,dtype=dtype,method=method)))
    end
    rewiring_fractions = append!(["0.0"],rewiring_fractions)

    return data,rewiring_fractions
end

function get_triangle_diffusions_plotting_data(gname::String,
                                    beta::Float64=0.05,
                                    gamma::Float64=0.05,
                                    dloc::String="pipeline/data/",
                                    dtype::String="tinfs",
                                    method::String="seir",
                                    ntrials::Int=50)

    A = loadGraph(gname,"pipeline/graphs/")
    nnodes = size(A,1) 

    data,rewiring_fractions = load_triangle_diffusion_data(gname,beta,gamma,
                                    dloc,dtype,method)
    
    if dtype == "tinfs"
        clabel = "Fraction of Infected Nodes"
        data = hcat(map(x->last.(x),data)...)./nnodes
    elseif dtype == "cinfs"
        clabel = "Normalized Maximum Infection Size"
        data = hcat(map(x->maximum.(x),data)...)./nnodes
    end

    nrows,ncols = size(data)
    nrows = Int(nrows/ntrials)
    aggreated_data = zeros(Float64,nrows,ncols)    
    for col in eachindex(1:ncols)
        for row in eachindex(1:nrows)
            offset = (row-1)*ntrials
            aggreated_data[row,col] = sum(data[offset+1:offset+ntrials,col])/ntrials
        end
    end
    return data,aggreated_data,clabel,rewiring_fractions
end


function get_rewired_data(gnames::Vector{String})
    result = Dict{String,Dict{Float64,Array}}()
    println("loading data")
    @showprogress for gname in gnames
        println("working on $gname")
        betas = get_betas(gname)
        gname_data = Dict{Float64,Array}()
        for beta in betas
            #load data and append 
            gname_data[beta] = get_plotting_data(gname,beta=beta)[2]
        end
        result[gname] = gname_data
    end
    return result
end

get_rewired_data(gname::String) = get_rewired_data([gname])


