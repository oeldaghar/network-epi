working_dir = "/p/mnt/scratch/network-epi/"
#set up 
include(joinpath(working_dir,"code","graph-io.jl"))

#sets up directories
println("setting things up") 
include(joinpath(working_dir,"input","setup.jl"))

#rewire base graphs 
println("Performing graph rewirings")
include(joinpath(working_dir,"code","rewiring.jl"))

### ncp computaions 
#acl ncp
println("working on acl ncp")
include(joinpath(working_dir,"code","ncp","ncp-acl-script.jl"))

#epidemic ncp
println("working on epidemic ncp")
include(joinpath(working_dir,"code","ncp","epidemic-ncp-script.jl"))

### diffusions
#uniform diffusions
println("working on uniform diffusions")
include(joinpath(working_dir,"code","fast-qdiffusions.jl"))

#triange weighted 
# println("working on triangle weighted diffusions")
# include("fast-weighted-diffusions.jl")

#generate all figures


