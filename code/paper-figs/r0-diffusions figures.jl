# r0-diffusions figures
using DelimitedFiles
using StatsBase
using Plots 
using Measures
mainDir = "/p/mnt/scratch/network-epi/"
include(joinpath(mainDir,"code","graph-io.jl"))
include(joinpath(mainDir,"code","data-io.jl"))

ENV["GKSwstype"] = "100"

gnames = [
    # "commutes-all","mexico", "filtered", #row1
    # "cn-moduillinois", "cn-Penn", "cn-modWiscon", #row2...
    # "dblp","enron","anon", #row 3
    # "cit-HepPh", "slashdot", "flickr", 
    # "geometric",
    # "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.9-0.0-100-5",
    "study-11-2023-0-noweak.smat",
    # "study-11-2023-1-longrange-1.smat",
    # "study-11-2023-1-longrange-2.smat",
    # "study-11-2023-1-longrange-3.smat",
    # "study-11-2023-1-longrange-5.smat",
    # "study-11-2023-1-longrange-8.smat",
    # "study-11-2023-1-longrange-10.smat",
    # "study-11-2023-1-longrange-12.smat",
    # "study-11-2023-1-longrange-15.smat",
    "study-11-2022-1.smat",
    "study-11-2022-10.smat",
    "study-11-2022-20.smat",
    "study-11-2022-30.smat",
    "study-11-2022-40.smat", 
    "study-11-2022-50.smat",
    # "covidflows-2020_08_31-filtered-20.smat",
    # "commutes-all.smat",
    # "modmexico-city.smat",
    # "email-Enron.smat",
]


gnames = getgnames("noweak","input/graphs/")
push!(gnames,getgnames("11-2022","tmp-study")...)

gnames = getgnames("study-25","input/graphs/")


# gnames = [gname]
# R0s = [1.1,2,5,10,20,50]#,100,200,500,1000]
# R0s = [1.1,2,5,10,20,25,30,35,40,45,50,100]#,200,500,1000]
# R0s = [25.0,30,35,40,45]
R0s = [1.1, 1.5, 2.5, 3, 3.5, 4, 4.5, 6, 7, 8, 9,
10,20,25,35,35,40,45,50,60]#, 70, 75, 80, 90,
# 100,200,500,1000] 

R0s = [1.1, 5, 10, 20, 25, 26, 27, 28, 29, 30, 35, 40, 45, 50, 60, 75]
eig_data = readdlm("pipeline/data/dominant-eigvals-0.txt")
eig_data = Dict(zip(eig_data[:,1],eig_data[:,2]))


function beta_from_r0(r0,lam,gamma::Float64=5e-2)
    #r0 = lam*beta/gamma
    return r0*gamma/lam
end

function load_diffusion_data(gname::String,eig_data::Dict;
          R0s::Vector{Float64}=[2.0],
          method::String="seir")

  # rps = ["10000.0", "1000.0", "100.0", "50.0", "10.0", "5.0",]
  gnames = getgnames(gname,"pipeline/graphs/")#["rewired-10000.0-$gname"; gname; "er-10000.0-$gname"]
  # gnames = [gname]
  g = canonical_graph_name(gname)

  #load in data   
  data = Dict{String,Dict}() #gname -> beta -> total
  dloc = "pipeline/data/$(g[1:end-5])/diffusions/uniform/scratch/"

  for h in gnames
    h_data = Dict{Float64,Vector{Float64}}() 
    #get betas 
    h_lambda = eig_data[h]
    betas = [beta_from_r0(x,h_lambda) for x in R0s]
    inds =  betas.<=1

    betas = betas[inds]
    R0s = R0s[inds]
    
    for (r0,beta) in zip(R0s,betas)
      h_data[r0] = sum.(read_inf_data(h,dloc=dloc,beta=beta,dtype="tinfs",method=method))
    end
    data[h] = deepcopy(h_data)
  end

  return data
end

# gnames = getgnames("longrange-5-","input/graphs/")
#testing 
load_diffusion_data(gnames[1],eig_data,R0s=[1.1])

#load in data
data = Dict{String,Dict}()
@showprogress for (ind,gname) in enumerate(gnames)
  tmp = load_diffusion_data(gname,eig_data,R0s = R0s)
  for k in keys(tmp)
    data[k] = tmp[k]
  end
end
#aggreate data 
aggregated_data = Dict{String,Dict}()
@showprogress for gname in gnames
  #load in graph once 
  nnodes = graphSize(gname)

  # graph_variants = ["rewired-10000.0-$gname"; gname; "er-10000.0-$gname"]
  graph_variants = getgnames(gname,"pipeline/graphs/")
  for h in graph_variants
    h_data = data[h]
    #process betas for a single graph 
    tmp = Dict()
    for r0 in keys(h_data)
      tmp[r0] = aggregate_diffusion_data(hcat(h_data[r0]))./nnodes
    end
    #add to target 
    aggregated_data[h] = deepcopy(tmp)
  end
end 

function qimpact_metric(data)
    return 1-( mean(data[2:end])/(data[1]) )
end
  
function make_r0_plot(gname,aggregated_data)
  g_r0s = sort(collect(keys(aggregated_data[gname])))
  ys = [qimpact_metric(aggregated_data[gname][r0]) for r0 in g_r0s]
  base_total_infs = [aggregated_data[gname][r0][1] for r0 in g_r0s]
  avg_q_total_infs = [mean(aggregated_data[gname][r0][2:end]) for r0 in g_r0s]

  f = plot(g_r0s,ys,label="Qimpact",
    xlims=(0.8,1.2*maximum(g_r0s)),xscale=:log10,
    ylims=(-0.03,1.03),
    ylabel="Node Proportion",
    markerstrokewidth=0,markershape=:circle,
    linestyle=:dash,)
  plot!(f,g_r0s,base_total_infs,color=:red,label="TotalInfs",
    xlabel="R0",markerstrokewidth=0,markershape=:circle,
    linestyle=:dash,
    legend=:topleft)
  plot!(f,g_r0s,avg_q_total_infs,color=:green,label="AvgQInfs",
    xlabel="R0",markerstrokewidth=0,markershape=:circle,
    linestyle=:dash,
    legend=:topleft)
  title!(f,gname[1:end-5])
  plot!(f,xticks=(R0s,string.(R0s)),xrotation=90)
end

figs = []

for gname in gnames 
  f = make_r0_plot("$gname",aggregated_data)
  push!(figs,f)
end

make_r0_plot("er-10000.0-$(gnames[1])",aggregated_data)

# for i=1:lastindex(gnames)
#   gname = gnames[i]
#   f = figs[i]
#   plot!(f,dpi=250)
#   Plots.savefig(f,"/p/mnt/scratch/network-epi/code/paper-figs/qimpact-figs/$(gname[1:end-5]).png")
# end
figs
figs[1]
figs[2]
figs[3]
figs[4]
figs[5]
figs[6]
figs[7]
figs[8]
figs[9]
figs[10]
figs[11]
figs[12]
figs[13]
figs[14]
figs[15]
figs[16]
figs[17]
figs[18]
figs[19]
figs[20]
figs[21]
figs[22]


beta_from_r0(70,eig_data[gnames[end]])
beta_from_r0(30,eig_data[gnames[2]])



figs[1]
figs[2]
figs[11]
figs[21]
figs[31]
figs[41]
figs[51]

f = make_r0_plot("er-10000.0-$gname",aggregated_data)

tmpfs = deepcopy(figs)
ylabel!.(tmpfs,"")
xlabel!.(tmpfs,"")
# title!.(tmpfs,"")
plot!.(tmpfs,legend=false)
newf = plot(tmpfs[[1;5;7;9;1;10;11;12]]...,layout=(2,4),
        margins=(-0.5Measures.mm),size=(1500,600))


gname = getgnames("filtered","input/graphs/")[1]
g_r0s = sort(collect(keys(aggregated_data[gname])))
ys = [qimpact_metric(aggregated_data[gname][r0]) for r0 in g_r0s]
base_total_infs = [aggregated_data[gname][r0][1] for r0 in g_r0s]


#making plots for PPR AUC versus qimpact for( networks 
function area_impact_plot(area_data,impact_data)
  xs = log10.(Float64.(area_data))
  ys = impact_data

  p = scatter(10 .^xs, ys,
        leg=false,
        color=:red3,
        xscale=:log10,
        markerstrokewidth=0,
        ylabel="Quarantine Impact")
  return p 
end

#from local structure error file 
# log_auc_measure
# load_ncps(gnames,ncptype="acl") function 
#plot against each other 

#ppr data 
gnames = getgnames("study-25","input/graphs/")

acl_ncpdata = load_ncps(gnames,ncptype="acl")
nnodes_dict = Dict()
@showprogress for gname in gnames
  @show(gname)
  A = loadGraph(gname,"input/graphs/")
  nnodes_dict[gname] = size(A,1)
end

ncp_area = log_auc_measure(acl_ncpdata,graphSizes=nnodes_dict)


ys = [qimpact_metric(aggregated_data[gname][35]) for gname in gnames]
xs = [ncp_area[gname] for gname in gnames]

scatter(xs,ys,leg=false)


f = area_impact_plot(xs,ys)
x,y = acl_ncpdata[gnames[1]]

#get mins 
minx,miny = get_approximate_mins(x,y)

#extrapolate
nnodes = nnodes_dict[canonical_graph_name(gnames[1])]
extrapolate!(minx,miny,nnodes)

#compute auc from y=1 to conductance minimum
auc = approximate_auc(log10.(minx), abs.(log10.(miny))) #auc from y=1 to min conductance
auc/(log10(nnodes/2))
gnames[1]


#for use with study-24 + study-25 graphs
function r0_qpercent_contours(gname::String,aggregated_data::Dict;
        r0s::Vector{Float64} = [1.1, 5, 10, 20, 25, 26, 27, 28, 29, 30, 35, 40, 45, 50, 60, 75],
        qpercents::Vector{Int}=collect(0:15))

    figs = []

    #main loop 
    @showprogress for r0 in r0s
      
      rps = ["10000.0","1000.0","100.0","50.0","25.0","10.0","5.0","1.0","0.5","0.1","0.05","0.01"]
      sorted_gnames = vcat(["rewired-$x-$(gname)" for x in rps],[gname], ["er-$x-$(gname)" for x in reverse(rps)])

      #aggreagted_data[gname] = Dict(r0->diffusion data)
      data = hcat([aggregated_data[g][r0] for g in sorted_gnames]...)
      f = Plots.contourf(1:25,1:16,(x,y)->data[y,x],
            c=:heat,
            clims=(0,1),
            levels=0.0:0.2:1.0,
            xrotation=90)
      
      Plots.plot!(f,xticks=(1:25, vcat(rps,[0],reverse(rps))),margins=9Measures.mm)
      
      push!(figs,f)
    end    
    return figs,r0s
end


### Portion for making figures for R0 diffusions - this is the spatial figure in the paper (GeometricCommunities) ###
gname = "study-25-150.smat"

gs = ["study-25-1.smat","study-25-2.smat","study-25-150.smat"]
figs,rs = r0_qpercent_contours(gs[3],aggregated_data)


f = figs[8]
plot!(f,colorbar=false,xticks=false,yticks=false,margins=0Measures.mm,
    dpi=500)

#make highresolution figures 
for gname in gs 
  figs,rs = r0_qpercent_contours(gname,aggregated_data)
  f = figs[8] #r0=28
  plot!(f,colorbar=false,xticks=false,yticks=false,margins=0Measures.mm,
    dpi=500)
  Plots.savefig(f,joinpath(mainDir,"code","paper-figs","spatial-figs","uplot-r0-28-$(gname[1:end-5])"))
end