#set up directories and generate rewired graphs.
#after this script, one needs to run the diffusions and get ncpdata
parentDir = "/p/mnt/scratch/network-epi/"
cd(parentDir)
graphs = readdir("input/graphs/")
filter!(x->endswith(x,".smat"),graphs)

#=
set up directories
/pipeline
	graphs/ #stores all graphs+rewired versions
	data/
		graphname/
			imgs/
				missed-sets/ #beta,qpercent,gname ~19*16*25 imgs per base graph
			ncpdata/
			diffusions/
				uniform/
					plotting-data/
					scratch/
				triangle-weighted/
					plotting-data/
				new-triangle-weighted/
					plotting-data/
=#

graphDir = joinpath(parentDir,"pipeline/graphs/")
if !ispath(graphDir)
	mkpath(graphDir)
end

for graph in graphs
	#set up directories
	ncpDir = joinpath(parentDir,"pipeline/data/$(graph[1:end-5])/ncpdata/")
	if !ispath(ncpDir)
		mkpath(ncpDir)
	end
	uniformDir = joinpath(parentDir,"pipeline/data/$(graph[1:end-5])/diffusions/uniform/plotting-data/")
	if !ispath(uniformDir)
		mkpath(uniformDir)
	end

	#for extra parameter testing. won't have data for other graphs 
	uniformDir = joinpath(parentDir,"pipeline/data/$(graph[1:end-5])/diffusions/uniform/scratch/")
	if !ispath(uniformDir)
		mkpath(uniformDir)
	end

	#for study + lfr graph testing 
	uniformDir = joinpath(parentDir,"pipeline/data/$(graph[1:end-5])/diffusions/uniform/new/")
	if !ispath(uniformDir)
		mkpath(uniformDir)
	end

	triangleDir = joinpath(parentDir,"pipeline/data/$(graph[1:end-5])/diffusions/triangle-weighted/plotting-data/")
	if !ispath(triangleDir)
		mkpath(triangleDir)
	end
	#for new triangle experiments. condense later.
	triangleDir = joinpath(parentDir,"pipeline/data/$(graph[1:end-5])/diffusions/new-triangle-weighted/plotting-data/")
	if !ispath(triangleDir)
		mkpath(triangleDir)
	end

	imgDir = joinpath(parentDir,"pipeline/data/$(graph[1:end-5])/imgs/missed-sets/")
	if !ispath(imgDir)
		mkpath(imgDir)
	end
	#copy original graph
	if !isfile(joinpath(parentDir,"pipeline/graphs/$graph"))
		cp(joinpath(parentDir,"input/graphs/",graph),joinpath(parentDir,"pipeline/graphs/$graph"))
	end
end


