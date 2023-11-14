#general functions for loading/writing data
using SparseArrays 
using MatrixNetworks
using LinearAlgebra
"""
    loadGraph(gname::String,gpath::String="";sym::Bool=true)

load in an SMAT representation of a graph, remove the diagonal and take the largest component.
use sym=true to symmetrize the graph before taking largest component.
"""
function loadGraph(gname::String,gpath::String="";sym::Bool=true)
    A = MatrixNetworks.readSMAT(gpath*gname)
    @assert(isequal(size(A)...))
    if sym
        A = max.(A,A')
    end
    A.-=spdiagm(0=>diag(A))
    dropzeros!(A)
    fill!(A.nzval,1.0)
    A = largest_component(A)[1]
    return A
end

"""
    getgnames(gname::String,gpath::String)

function for getting graphs names based on case-insenitive substring
graph names are sorted by legnth so that the base graph occurs first 
"""
function getgnames(gname::String,gpath::String)
    return sort(filter(x->occursin(lowercase(gname),lowercase(x)) && endswith(x,".smat"),readdir(gpath)),by=length)
end

"""
    include_graph(gname,gpath)

get rewired graph names from input graph name but exclude sparse graphs
"""
function include_graph(gname,gpath)
    gnames = getgnames(gname,gpath)
    if !startswith(gname,"cn-") #CA-Astro picks up cn-CA-Astro. want to filter these out.
        filter!(x->!occursin("cn-",x),gnames)
    end
    return gnames 
end

"""
    get_graph_stats(fname::String;gpath::String="")

peek at number of nodes and edges a graph has before loading it in
"""
function get_graph_stats(fname::String;gpath::String="input/graphs/")
    open(joinpath(gpath,fname),"r") do io
        nnodes,nedges = parse.(Int,split(readline(io),(' ','\t')))[2:3]
        return nnodes,nedges
    end
end


"""
    writeSMAT(A::SparseMatrixCSC,fname::String)

store A using an SMAT representation
"""
function writeSMAT(A::SparseMatrixCSC,fname::String)
    Is,Js,Vs = findnz(A)
    nedges = length(Is)
    open(fname,"w") do file
        write(file,"$(size(A,1)) $(size(A,2)) $nedges\n")
        for i=1:nedges
            write(file,"$(Is[i]-1) $(Js[i]-1) $(Vs[i])\n") #apparently smat assumes node labels start at 0
        end
    end
end

"""
    canonical_graph_name(gname::String)

return the name of the base graph used to generate gname
"""
function canonical_graph_name(gname::String)
    if startswith(gname,"er-") || startswith(gname,"rewired-")
        g = join(split(gname,"-")[3:end],"-")
    elseif startswith(gname,"triangle-rewired")
        g = join(split(gname,"-")[5:end],"-")
    else
        g = gname
    end
    return g
end

#path specific 
function graphSize(gname)
    A = loadGraph(joinpath("pipeline/graphs/",gname))
    return size(A,1)
end
  