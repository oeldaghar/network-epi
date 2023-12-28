#generate a table of graph statistics for all graphs.
#information should include nnodes,nedges,avgd,deepest kcore, %nodes in 2-core, spectral radius

using LinearAlgebra
using MatrixNetworks
using Arpack
using DelimitedFiles

#allows for execution from command line as well as an ide so files can be run modularlly
mainDir = joinpath("/",split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))]...)

include(joinpath(mainDir,"code/graph-io.jl"))

#only do this for base graphs 
gpath = "input/graphs/"

gnames = filter(x->endswith(x,".smat"),readdir(gpath))

function get_graph_info(gnames::Vector{String},corenum::Int=3)
    #load graph, compute stats, and return output 
    data = Vector{Vector}()
    header = vec(["gname" "nnodes" "nedges" "avgd" "max kcore" "small_core_frac" "spectral radius"])
    push!(data,header)
    
    #load precomputed eigenvalue data
    eig_data = readdlm("pipeline/data/dominant-eigvals-0.txt")

    for gname in gnames 
        A = loadGraph(gname,"pipeline/graphs/")
        @assert(issymmetric(A))

        #size info
        nnodes = size(A,1)
        nedges = nnz(A)
        avgd = round(nedges/nnodes,digits=1)

        #kcore info
        kcore = corenums(A)[1]
        max_kcore = maximum(kcore)

        small_core_frac = round(sum(kcore.<corenum)/nnodes,digits=2)

        #spectral radius 
        max_eigval = round(eig_data[findfirst(eig_data[:,1].==gname),2],digits=2)
        
        
        push!(data,vec([gname[1:end-5] nnodes nedges avgd max_kcore small_core_frac max_eigval]))
    end
    return data 
end

data = get_graph_info(["study-25-1.smat","study-25-2.smat","study-25-150.smat"])

#format as latex table 
#\begin{table}
# \centering
# \begin{tabular}{c|c}
#      &  \\
#      & 
# \end{tabular}
# \caption{Caption}
# \label{tab:my_label}
# \end{table}

#will do formatting for the tabular part only. so use "&" as delimiter, and end lines with \\
# format each row and append to file. final edits will be made in latex

fname = joinpath(mainDir,"code","paper-figs","graph-stats","graph-statistics-table.txt")
open(fname,"w") do io
    for row in data 
        newrow = join(row," & ")*" \\\\ "
        writedlm(io, [newrow])
    end
end

