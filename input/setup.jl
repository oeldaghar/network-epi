#set up directories and generate rewired graphs.
#after this script, one needs to run the diffusions and get ncpdata
mainDir = joinpath("/",split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))]...)
cd(mainDir)
graphs = readdir(joinpath(mainDir,"input/graphs/"))
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

graphDir = joinpath(mainDir,"pipeline/graphs/")
if !ispath(graphDir)
	mkpath(graphDir)
end

for graph in graphs
	#set up directories
	ncpDir = joinpath(mainDir,"pipeline/data/$(graph[1:end-5])/ncpdata/")
	if !ispath(ncpDir)
		mkpath(ncpDir)
	end
	uniformDir = joinpath(mainDir,"pipeline/data/$(graph[1:end-5])/diffusions/uniform/plotting-data/")
	if !ispath(uniformDir)
		mkpath(uniformDir)
	end

	#for extra parameter testing. won't have data for other graphs 
	uniformDir = joinpath(mainDir,"pipeline/data/$(graph[1:end-5])/diffusions/uniform/scratch/")
	if !ispath(uniformDir)
		mkpath(uniformDir)
	end

	#for study + lfr graph testing 
	uniformDir = joinpath(mainDir,"pipeline/data/$(graph[1:end-5])/diffusions/uniform/new/")
	if !ispath(uniformDir)
		mkpath(uniformDir)
	end

	triangleDir = joinpath(mainDir,"pipeline/data/$(graph[1:end-5])/diffusions/triangle-weighted/plotting-data/")
	if !ispath(triangleDir)
		mkpath(triangleDir)
	end
	#for new triangle experiments. condense later.
	triangleDir = joinpath(mainDir,"pipeline/data/$(graph[1:end-5])/diffusions/new-triangle-weighted/plotting-data/")
	if !ispath(triangleDir)
		mkpath(triangleDir)
	end

	imgDir = joinpath(mainDir,"pipeline/data/$(graph[1:end-5])/imgs/missed-sets/")
	if !ispath(imgDir)
		mkpath(imgDir)
	end
	#copy original graph
	if !isfile(joinpath(mainDir,"pipeline/graphs/$graph"))
		cp(joinpath(mainDir,"input/graphs/",graph),joinpath(mainDir,"pipeline/graphs/$graph"))
	end
end


