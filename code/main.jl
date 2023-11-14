#allows for execution from command line as well as an ide so files can be run modularlly
mainDir = joinpath(split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))])

#set up 
include(joinpath(mainDir,"code","graph-io.jl"))

#sets up directories
println("setting things up") 
include(joinpath(mainDir,"input","setup.jl"))

#rewire base graphs 
println("Performing graph rewirings")
include(joinpath(mainDir,"code","rewiring.jl"))

### ncp computaions 
#acl ncp
println("working on acl ncp")
include(joinpath(mainDir,"code","ncp","ncp-acl-script.jl"))

#epidemic ncp
println("working on epidemic ncp")
include(joinpath(mainDir,"code","ncp","epidemic-ncp-script.jl"))

### diffusions
#uniform diffusions
println("working on uniform diffusions")
include(joinpath(mainDir,"code","fast-qdiffusions.jl"))

#triange weighted 
# println("working on triangle weighted diffusions")
# include("fast-weighted-diffusions.jl")

#generate all figures


