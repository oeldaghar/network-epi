parent_dir = "/p/mnt/scratch/network-epi/"
include(joinpath(parent_dir,"code/graph-io.jl"))
include(joinpath(parent_dir,"code/data-io.jl"))
include(joinpath(parent_dir,"code/fast-diffusion.jl"))
# include(joinpath(parent_dir,"code/rt.jl"))

using SparseArrays
using DelimitedFiles
using Plots
using Measures
ENV["GKSwstype"] = "100"

#(dominant RAW/PL eig ratio) vs (Qpercent for total infs <20%)

#for all graphs, compute and save eigenvalue information
using Arpack

function get_dominant_eigenvalue(A::SparseMatrixCSC)
    lam,v = eigs(A,nev=1,maxiter=1000)
    lam = lam[1]
    return lam,norm(A*v-lam*v)
end

function compute_eigenvalues(gnames::Vector{String};rseed::Int=0)
    if isfile("pipeline/data/dominant-eigvals-$rseed.txt")
        data = readdlm("pipeline/data/dominant-eigvals-$rseed.txt")
    else
        data = Vector()
    end

    @showprogress for gname in gnames
        #check if we already computed this value
        if !(gname in data[:,1])
            Random.seed!(rseed)
            println("working on lambda1 for $gname")
            A = loadGraph(gname,"pipeline/graphs/")
            lam,err = get_dominant_eigenvalue(A)
            open("pipeline/data/dominant-eigvals-$rseed.txt","a") do io
                writedlm(io,[gname lam err])
            end
        end
    end
    return true
end

gnames = readdir(joinpath(parent_dir,"pipeline/graphs/"))
filter!(x->endswith(x,".smat"),gnames)
compute_eigenvalues(gnames)

#lam1(pl)/lam1(orig)
function get_lambda_plotting_data(gnames::Vector{String})
    lams = readdlm("pipeline/data/dominant-eigvals-0.txt")

    xs = Vector{Vector{Float64}}()
    ys = Vector{Vector{Float64}}()
    @showprogress for gname in gnames
        println("working on graph $gname")
        xbeta = Vector{Float64}()
        ybeta = Vector{Float64}()
        for beta in vcat(1e-3:1e-3:1e-2,2e-2:1e-2:1e-1)
            aggregated_data = get_plotting_data(gname,beta=beta)[2]
            
            #get ratio of total infs
            ncols = size(aggregated_data)[2]
            orig_col = ceil(Int,ncols/2)
            push!(xbeta,sum(aggregated_data[2:end,1])/sum(aggregated_data[2:end,orig_col]))
            
            #get ratio of dominant eigvals
            raw_lam = lams[findfirst(lams[:,1].==gname),2]
            pl_lam = lams[findfirst(lams[:,1].=="rewired-10000.0-$(gname)"),2]

            push!(ybeta,pl_lam/raw_lam)
        end
        push!(xs,xbeta)
        push!(ys,ybeta)
    end

    xs = hcat(xs...)
    ys = hcat(ys...)
    return xs,ys    
end

function get_lambda_plotting_data1(gnames::Vector{String})
    lams = readdlm("pipeline/data/dominant-eigvals-0.txt")

    xs = Vector{Vector{Float64}}()
    ys = Vector{Vector{Float64}}()
    @showprogress for gname in gnames
        #original graph
        A = loadGraph(gname,"input/graphs/")
        orig_deg = vec(sum(A;dims=2))
        orig_dmax = maximum(orig_deg)
        orig_davg = mean(orig_deg)

        #PL rewired graph 
        A = loadGraph("rewired-10000.0-$gname","pipeline/graphs/")
        pl_deg = vec(sum(A;dims=2))
        pl_dmax = maximum(pl_deg)
        pl_davg = mean(pl_deg)

        println("working on graph $gname")
        xbeta = Vector{Float64}()
        ybeta = Vector{Float64}()
        for beta in vcat(1e-3:1e-3:1e-2,2e-2:1e-2:1e-1)
            aggregated_data = get_plotting_data(gname,beta=beta)[2]
            
            #get ratio of total infs
            ncols = size(aggregated_data)[2]
            orig_col = ceil(Int,ncols/2)
            push!(xbeta,sum(aggregated_data[2:end,1])/sum(aggregated_data[2:end,orig_col]))
            
            #get ratio of dominant eigvals
            raw_lam = lams[findfirst(lams[:,1].==gname),2]
            pl_lam = lams[findfirst(lams[:,1].=="rewired-10000.0-$(gname)"),2]

            raw_statistic = (raw_lam-sqrt(orig_dmax))/(orig_dmax-sqrt(orig_dmax))
            pl_statistic = (pl_lam-sqrt(pl_dmax))/(pl_dmax-sqrt(pl_dmax))
            push!(ybeta,pl_statistic/raw_statistic)
        end
        push!(xs,xbeta)
        push!(ys,ybeta)
    end

    xs = hcat(xs...)
    ys = hcat(ys...)
    return xs,ys    
end

function get_lambda_plotting_data2(gnames::Vector{String})
    lams = readdlm("pipeline/data/dominant-eigvals-0.txt")

    xs = Vector{Vector{Float64}}()
    ys = Vector{Vector{Float64}}()
    @showprogress for gname in gnames
        #original graph
        A = loadGraph(gname,"input/graphs/")
        orig_deg = vec(sum(A;dims=2))
        orig_dmax = maximum(orig_deg)
        orig_davg = mean(orig_deg)

        #PL rewired graph 
        A = loadGraph("rewired-10000.0-$gname","pipeline/graphs/")
        pl_deg = vec(sum(A;dims=2))
        pl_dmax = maximum(pl_deg)
        pl_davg = mean(pl_deg)

        println("working on graph $gname")
        xbeta = Vector{Float64}()
        ybeta = Vector{Float64}()
        for beta in vcat(1e-3:1e-3:1e-2,2e-2:1e-2:1e-1)
            aggregated_data = get_plotting_data(gname,beta=beta)[2]
            
            #get ratio of total infs
            ncols = size(aggregated_data)[2]
            orig_col = ceil(Int,ncols/2)
            
            pl_quarantine_impact = sum(aggregated_data[2:end,1])/aggregated_data[1,1]
            orig_quarantine_impact = sum(aggregated_data[2:end,orig_col])/aggregated_data[1,orig_col]
            
            push!(xbeta,pl_quarantine_impact/orig_quarantine_impact)
            
            #get ratio of dominant eigvals
            raw_lam = lams[findfirst(lams[:,1].==gname),2]
            pl_lam = lams[findfirst(lams[:,1].=="rewired-10000.0-$(gname)"),2]

            # raw_statistic = (raw_lam-sqrt(orig_dmax))/(orig_dmax-sqrt(orig_dmax))
            # pl_statistic = (pl_lam-sqrt(pl_dmax))/(pl_dmax-sqrt(pl_dmax))
            raw_statistic = raw_lam 
            pl_statistic = pl_lam

            push!(ybeta,pl_statistic/raw_statistic)
        end
        push!(xs,xbeta)
        push!(ys,ybeta)
    end

    xs = hcat(xs...)
    ys = hcat(ys...)
    return xs,ys    
end

#first plot 
#PL/RAW extremal eigvals for all graphs 


function make_eigenvalue_plot(sname::String="study")
    
    gnames = ["study-11-2022-5.smat","study-11-2022-10.smat",
        "study-11-2022-15.smat","study-11-2022-20.smat","study-11-2022-25.smat",
        "study-11-2022-30.smat","study-11-2022-35.smat","study-11-2022-40.smat",
        "study-11-2022-45.smat","study-11-2022-50.smat"]

    xs,ys = get_lambda_plotting_data(gnames)

    f = scatter(xs[:,1],ys[:,1],leg=false,
            xlims=(0.5,maximum(35,maximum(xs)+5)),ylims=(0.53,0.68),
            # alpha=(11)/20,
            xscale=:log10,
            markerstrokewidth=0,
            xlabel="Total Infs Ratio (PL/Orig)",
            ylabel="Dominant Eigenvalue Ratio (PL/Orig)")
    for ind = 2:10
        scatter!(f,xs[:,ind],ys[:,ind],leg=false,
            # alpha=(10+ind)/20,
            markerstrokewidth=0)
    end

    plot!(f,[1;1],[ylims()...])
    plot!(f,[1;100],[1;1],fillrange=[0; 0],
            alpha=0.15,color=3)

    annotate!(f,0.8,0.675,text("5", 8))
    annotate!(f,0.8,0.655,text("10", 8))
    annotate!(f,0.95,0.643,text("15", 8))
    annotate!(f,0.9,0.625,text("20", 8))
    annotate!(f,0.9,0.615,text("25", 8))
    annotate!(f,0.9,0.568,text("30", 8))
    annotate!(f,0.9,0.558,text("35", 8))
    annotate!(f,0.88,0.55,text("40", 8))
    annotate!(f,22,0.575,text("45", 8))
    annotate!(f,18,0.582,text("50", 8))
    title!(f,"SEIR(β,0.05) - Study-11-2022-XXXX")

    annotate!(f,1.6,0.6,text("⟶ Larger Eigval\n     Less Infs for Original", 9))
    return f
end

f = make_eigenvalue_plot()
plot!(f,xlims=(0.5,220))
f

gnames = ["study-11-2022-5.smat","study-11-2022-10.smat",
"study-11-2022-15.smat","study-11-2022-20.smat","study-11-2022-25.smat",
"study-11-2022-30.smat","study-11-2022-35.smat","study-11-2022-40.smat",
"study-11-2022-45.smat","study-11-2022-50.smat"]

xs,ys = get_lambda_plotting_data(gnames)

extrema(xs)
extrema(ys)

f = scatter(xs[:,1],ys[:,1],leg=false,
    xlims=(0.5,max(35,1.2*maximum(xs)+5)),ylims=(0.53,0.68),
    # alpha=(11)/20,
    xscale=:log10,
    markerstrokewidth=0,
    xlabel="Total Infs Ratio (PL/Orig)",
    ylabel="Dominant Eigenvalue Ratio (PL/Orig)")
for ind = 2:10
    scatter!(f,xs[:,ind],ys[:,ind],leg=false,
    # alpha=(10+ind)/20,
    markerstrokewidth=0)
end

plot!(f,[1;1],[ylims()...])
plot!(f,[1;xlims(f)[2]],[1;1],fillrange=[0; 0],
    alpha=0.15,color=3)

annotate!(f,0.8,0.675,text("5", 8))
annotate!(f,0.8,0.655,text("10", 8))
annotate!(f,0.95,0.643,text("15", 8))
annotate!(f,0.9,0.625,text("20", 8))
annotate!(f,0.85,0.615,text("25", 8))
annotate!(f,0.85,0.568,text("30", 8))
annotate!(f,0.85,0.558,text("35", 8))
annotate!(f,0.8,0.55,text("40", 8))
annotate!(f,4.2,0.575,text("45", 8))
annotate!(f,3,0.582,text("50", 8))
title!(f,"SEIR(β,0.05) - Study-11-2022-XXXX")

annotate!(f,2,0.6,text("⟶ Larger Eigval\n     Less Infs for Original", 9))


#different graphs 
gs = ["dblp","anon","mexico","brightkite","gowalla","commutes",
    "study-11-2022-25","epinions1","flickr","geometric",
    "slashdot","filtered","enron"]
gnames = vcat(getgnames.(gs,"input/graphs/")...)

xs,ys = get_lambda_plotting_data(gnames)

extrema(xs)
extrema(ys)

f = scatter(xs[:,1],ys[:,1],leg=false,
    xlims=(0.4,max(35,1.2*maximum(xs)+5)),
    ylims=(0.1,1.0),
    # alpha=(11)/20,
    xscale=:log10,
    markerstrokewidth=0,
    xlabel="Total Infs Ratio (PL/Orig)",
    ylabel="Dominant Eigenvalue Ratio (PL/Orig)")
for ind = 2:13
    scatter!(f,xs[:,ind],ys[:,ind],leg=false,
    # alpha=(10+ind)/20,
    markerstrokewidth=0)
end
f
plot!(f,[1;1],[ylims()...])
plot!(f,[1;xlims(f)[2]],[1;1],fillrange=[0; 0],
    alpha=0.15,color=3)

annotate!(f,40,ys[1,1],text(gs[1], 8))
annotate!(f,35,ys[1,2],text(gs[2], 8))
annotate!(f,200,ys[1,3],text(gs[3], 8))
annotate!(f,6,ys[1,4],text(gs[4], 8))
annotate!(f,4,ys[1,5],text(gs[5], 8))
annotate!(f,25,ys[1,6]-0.01,text(gs[6], 8))
annotate!(f,200,ys[1,7]+0.025,text(gs[7], 8))
annotate!(f,3.5,ys[1,8]-0.005,text(gs[8], 8))
annotate!(f,4,ys[1,9],text(gs[9], 8))
annotate!(f,300,ys[1,10]+0.025,text(gs[10], 8))
annotate!(f,2,ys[1,11],text(gs[11], 8))
annotate!(f,130,ys[1,12],text(gs[12], 8))
annotate!(f,0.7,ys[1,13],text(gs[13], 8))
annotate!(f,2.2,0.4,text("⟶ Larger Eigval\n     Less Infs for Original", 9))




#another plot 
#zoom in on a single graph
gnames = getgnames("geo","input/graphs/")[1:1] 
xs,ys = get_lambda_plotting_data2(gnames)

f = scatter([1],[xs[1,1]],leg=false,color=1,
    xlabel="Beta",
    ylabel="Total Infs Ratio (PL/Orig)")
for ind=2:lastindex(xs)
    scatter!(f,[ind],[xs[ind,1]],color=1)
end 
plot!(f,xticks=(1:19, vcat(1e-3:1e-3:9e-3,1e-2:1e-2:1e-1)),
    xrotation=90,
    title="$(gnames[1][1:end-5])")




gs = ["dblp","anon","mexico","brightkite","gowalla","commutes",
    "study-11-2022-25","epinions1","geometric",
    "slashdot","filtered","enron"]
gnames = vcat(getgnames.(gs,"input/graphs/")...)


xs,ys = get_lambda_plotting_data1(gnames)

extrema(xs)
extrema(ys)

f = scatter(xs[:,1],ys[:,1],leg=false,
    xlims=(0.4,max(35,1.2*maximum(xs)+5)),
    ylims=(0.0,1.1),
    # alpha=(11)/20,
    xscale=:log10,
    markerstrokewidth=0,
    xlabel="Total Infs Ratio (PL/Orig)",
    ylabel="Eigenvalue Ratio (PL/Orig)")
for ind = 2:12
    scatter!(f,xs[:,ind],ys[:,ind],leg=false,
    # alpha=(10+ind)/20,
    markerstrokewidth=0)
end
f
plot!(f,[1;1],[ylims()...])
plot!(f,[1;xlims(f)[2]],[1;1],fillrange=[0; 0],
    alpha=0.15,color=3)

annotate!(f,40,ys[1,1],text(gs[1], 8))
annotate!(f,35,ys[1,2],text(gs[2], 8))
annotate!(f,200,ys[1,3],text(gs[3], 8))
annotate!(f,0.8,ys[1,4]-0.03,text(gs[4], 8))
annotate!(f,4,ys[1,5],text(gs[5], 8))
annotate!(f,28,ys[1,6],text(gs[6], 8))
annotate!(f,200,ys[1,7]+0.03,text(gs[7], 8))
annotate!(f,3.5,ys[1,8]-0.005,text(gs[8], 8))
annotate!(f,25,ys[1,9]+0.04,text(gs[9], 8))
annotate!(f,2.,ys[1,10],text(gs[10], 8))
annotate!(f,125,ys[1,11],text(gs[11], 8))
annotate!(f,2.5,ys[1,12],text(gs[12], 8))
# annotate!(f,0.7,ys[1,13],text(gs[13], 8))
# annotate!(f,2.2,0.4,text("⟶ Larger Eigval\n     Less Infs for Original", 9))





# aggregated_data = get_plotting_data(gnames[1],beta=1e-3)[2]


gs = ["dblp","anon","brightkite","gowalla",
    "epinions1",
    "slashdot","enron"]
gnames = vcat(getgnames.(gs,"input/graphs/")...)

xs,ys = get_lambda_plotting_data2(gnames)

extrema(xs)
extrema(ys)
xs
f = scatter(xs[:,1],ys[:,1],leg=false,
    xlims=(min(1e-3,minimum(xs)),max(35,1.2*maximum(xs)+5)),
    ylims=(0.0,max(1.1,maximum(ys)+1e-1)),
    # alpha=(11)/20,
    xscale=:log10,
    markerstrokewidth=0,
    xlabel="Total Infs Ratio (PL/Orig)",
    ylabel="Eigenvalue Ratio (PL/Orig)")
for ind = 2:7
    scatter!(f,xs[:,ind],ys[:,ind],leg=false,
    # alpha=(10+ind)/20,
    markerstrokewidth=0)
end
f
plot!(f,[1;1],[ylims()...])
plot!(f,[1;xlims(f)[2]],[1;1],fillrange=[0; 0],
    alpha=0.15,color=3)


annotate!(f,40,ys[1,1],text(gs[1], 8))
annotate!(f,35,ys[1,2],text(gs[2], 8))
annotate!(f,200,ys[1,3],text(gs[3], 8))
annotate!(f,0.8,ys[1,4]-0.03,text(gs[4], 8))
annotate!(f,4,ys[1,5],text(gs[5], 8))
annotate!(f,28,ys[1,6],text(gs[6], 8))
annotate!(f,200,ys[1,7]+0.03,text(gs[7], 8))
annotate!(f,3.5,ys[1,8]-0.005,text(gs[8], 8))
annotate!(f,25,ys[1,9]+0.04,text(gs[9], 8))
annotate!(f,2.,ys[1,10],text(gs[10], 8))
annotate!(f,125,ys[1,11],text(gs[11], 8))
annotate!(f,2.5,ys[1,12],text(gs[12], 8))


gnames

#make plots based on this data 
#metrics for graph 
# - lam1 (dominant eigenvalue) 
# - lam1  / sqrt(dmax)
# - (lam1-sqrt(dmax))/(dmax- sqrt(dmax)) #normalize values to be in [0,1]

#metrics for impact of quarantining 
# - sum(quarantine infs)
# - sum(quarantine infs) / (nonquarantine infs) 
# - 1+norm(qinfs.-nonqinfs) #close to 1 means no change, larger means bigger change 


#iterate over data and generate whatever xs,ys we need 
gs = ["dblp","anon","mexico","brightkite","gowalla","commutes",
    "study-11-2022-25","epinions1","flickr","geometric",
    "slashdot","filtered","enron","cit-h","cit-pat"]
gnames = vcat(getgnames.(gs,"input/graphs/")...)
# gnames = getgnames("geo","input/graphs/")
data = get_rewired_data(gnames)

#
function qimpact_metric(inf_data::Matrix{Float64})
    raw_col_ind = ceil(Int,size(inf_data)[2]/2)
    pl_data = inf_data[:,1]
    raw_data = inf_data[:,raw_col_ind]

    # pl_stat = sum(pl_data[2:end])/pl_data[1]
    # raw_stat = sum(raw_data[2:end])/raw_data[1]

    #aggregating over quarantine information
    pl_stat = sum(pl_data[2:end])
    raw_stat = sum(raw_data[2:end])

    return pl_stat/raw_stat
end

function spectral_metric(pl_lam::Float64,pl_deg::Vector,raw_lam::Float64,raw_deg::Vector)
    pl_stat = pl_lam/mean(pl_deg)
    raw_stat = raw_lam/mean(raw_deg)
    return pl_lam/raw_lam
end


# dominant eigenvalues
xs = Vector{Vector}()
ys = Vector{Vector}()
lams = readdlm("pipeline/data/dominant-eigvals-0.txt")

#look at non-spatial graphs 
tmp_gs = ["dblp","anon","brightkite","gowalla",
"epinions1","flickr",
"slashdot","enron","cit-HepPh","cit-Patents"]
tmp_gnames = vcat(getgnames.(tmp_gs,"input/graphs/")...)
for gname in tmp_gnames
    println("working on graph $gname")
    x_gname= Vector{Float64}()
    y_gname= Vector{Float64}()

    #load graph and get spectral information 
    A = loadGraph(gname,"input/graphs/")
    raw_deg = vec(sum(A;dims=2))
    raw_lam = lams[findfirst(lams[:,1].==gname),2]

    #PL rewired graph 
    A = loadGraph("rewired-10000.0-$gname","pipeline/graphs/")
    pl_deg = vec(sum(A;dims=2))
    pl_lam = lams[findfirst(lams[:,1].=="rewired-10000.0-$(gname)"),2]  
    
    for beta in vcat(1e-3:1e-3:9e-3,1e-2:1e-2:1e-1)#keys(data[gname])
        #single graph single beta get statistic
        agg_data = data[gname][beta]c

        push!(x_gname,qimpact_metric(agg_data))
        push!(y_gname,spectral_metric(pl_lam,pl_deg,raw_lam,raw_deg))
    end
    push!(xs,x_gname)
    push!(ys,y_gname)
end

xs = hcat(xs...)
ys = hcat(ys...)

f = scatter(xs,ys,leg=false,xscale=:log10,markerstrokewidth=0)
# for (i,g) in enumerate(tmp_gnames)#enumerate(gs)
#     annotate!(f,20,ys[1,i],text(tmp_gs[i],5))
# end

annotate!(f,10,ys[1,1]+0.02,text(tmp_gs[1],10))
annotate!(f,20,ys[1,2]+0.02,text(tmp_gs[2],10))
annotate!(f,1,ys[1,3]-0.02,text(tmp_gs[3],10))
annotate!(f,0.5,ys[1,4],text(tmp_gs[4],10))
annotate!(f,0.65,ys[1,5]+0.01,text(tmp_gs[5],10))
annotate!(f,4,ys[1,6],text(tmp_gs[6],10))
annotate!(f,2,ys[1,7],text(tmp_gs[7],10))
annotate!(f,2,ys[1,8]-0.02,text(tmp_gs[8],10))
annotate!(f,4,ys[1,9]-0.02,text(tmp_gs[9],10))
annotate!(f,100,ys[1,10]+0.02,text(tmp_gs[10],10))

plot!(f,xlims=(3e-1,550),ylims=(0.2,0.95))
plot!(f,[1;1],[ylims()...])
plot!(f,[1;xlims(f)[2]],[1;1],fillrange=[0; 0],
    alpha=0.15,color=3)

xlabel!(f,"Total Quarantine Infs Ratio (PL/Raw) ")
ylabel!(f,"Spectral Ratio (PL/Raw) ")


#spatial graphs ----------------------------------
# dominant eigenvalues
xs = Vector{Vector}()
ys = Vector{Vector}()
lams = readdlm("pipeline/data/dominant-eigvals-0.txt")

#look at non-spatial graphs 
tmp_gs = ["mexico","commutes","study-11-2022-25","geometric","filtered"]

tmp_gnames = vcat(getgnames.(tmp_gs,"input/graphs/")...)
for gname in tmp_gnames
    println("working on graph $gname")
    x_gname= Vector{Float64}()
    y_gname= Vector{Float64}()

    #load graph and get spectral information 
    A = loadGraph(gname,"input/graphs/")
    raw_deg = vec(sum(A;dims=2))
    raw_lam = lams[findfirst(lams[:,1].==gname),2]

    #PL rewired graph 
    A = loadGraph("rewired-10000.0-$gname","pipeline/graphs/")
    pl_deg = vec(sum(A;dims=2))
    pl_lam = lams[findfirst(lams[:,1].=="rewired-10000.0-$(gname)"),2]  
    
    for beta in vcat(1e-3:1e-3:9e-3,1e-2:1e-2:1e-1)#keys(data[gname])
        #single graph single beta get statistic
        agg_data = data[gname][beta]

        push!(x_gname,qimpact_metric(agg_data))
        push!(y_gname,spectral_metric(pl_lam,pl_deg,raw_lam,raw_deg))
    end
    push!(xs,x_gname)
    push!(ys,y_gname)
end


xs = hcat(xs...)
ys = hcat(ys...)

f1 = scatter(xs,ys,leg=false,xscale=:log10,markerstrokewidth=0)
# for (i,g) in enumerate(tmp_gnames)#enumerate(gs)
#     annotate!(f,20,ys[1,i],text(tmp_gs[i],5))
# end

annotate!(f1,4,ys[1,1]+0.01,text(tmp_gs[1],10))
annotate!(f1,20,ys[1,2]+0.02,text(tmp_gs[2],10))
annotate!(f1,150,ys[1,3]-0.02,text(tmp_gs[3],10))
annotate!(f1,200,ys[1,4]-0.02,text(tmp_gs[4],10))
annotate!(f1,150,ys[1,5],text(tmp_gs[5],10))
plot!(f1,xlims=(3e-1,550),ylims=(0.53,0.87))
plot!(f1,[1;1],[ylims()...])


plot!(f1,[1;xlims(f1)[2]],[ylims(f1)[2];ylims(f1)[2]],fillrange=[0; 0],
    alpha=0.15,color=3)

xlabel!(f1,"Total Quarantine Infs Ratio (PL/Raw) ")
ylabel!(f1,"Spectral Ratio (PL/Raw) ")


plot!(f,dpi=500)
Plots.savefig(f,"code/paper-figs/eigenvalue-plots/eigval-plot.png")


plot!(f1,dpi=500)
Plots.savefig(f1,"code/paper-figs/eigenvalue-plots/eigval-plot-high-deviation.png")



ff = deepcopy(f)
ff1 = deepcopy(f1)
plot!(ff,xticks=[],xlabel="")
fig = Plots.plot(ff,ff1,layout=(2,1),size=(640,800),dpi=500,
        leftmargin=3Measures.mm)


Plots.savefig(fig,"code/paper-figs/eigenvalue-plots/eigval-plots-combined.png")


#look at individual graphs as we vary beta 
tmp_gname = gnames[5]
x_betas = Vector{Float64}()
for beta in keys(data[tmp_gname])
    push!(x_betas,qimpact_metric(data[tmp_gname][beta]))
end
betas = collect(keys(data[tmp_gname]))
beta_perm = sortperm(betas)

f = scatter(1:lastindex(x_betas),x_betas[beta_perm],leg=false,yscale=:log10)
plot!(f,xticks=(1:lastindex(x_betas), betas[beta_perm]),
    xrotation=90,
    title="$(tmp_gname[1:end-5])")

# qimpact_metric(data[gnames[2]][0.004])
lams[]
#make image for spectral gap 
using Arpack
A = loadGraph(getgnames("flickr","input/graphs/")[1],"input/graphs/")
lams,vs = eigs(A,nev=3,maxiter=1000)

lams


# ----------------------------------------------------------
# study plot
gnames = ["study-11-2022-5.smat","study-11-2022-10.smat",
"study-11-2022-15.smat","study-11-2022-20.smat","study-11-2022-25.smat",
"study-11-2022-30.smat","study-11-2022-35.smat","study-11-2022-40.smat",
"study-11-2022-45.smat","study-11-2022-50.smat"]

data = get_rewired_data(gnames)

xs = Vector{Vector}()
ys = Vector{Vector}()
lams = readdlm("pipeline/data/dominant-eigvals-0.txt")

#look at non-spatial graphs 
for gname in gnames
    println("working on graph $gname")
    x_gname= Vector{Float64}()
    y_gname= Vector{Float64}()

    #load graph and get spectral information 
    A = loadGraph(gname,"input/graphs/")
    raw_deg = vec(sum(A;dims=2))
    raw_lam = lams[findfirst(lams[:,1].==gname),2]

    #PL rewired graph 
    A = loadGraph("rewired-10000.0-$gname","pipeline/graphs/")
    pl_deg = vec(sum(A;dims=2))
    pl_lam = lams[findfirst(lams[:,1].=="rewired-10000.0-$(gname)"),2]  
    
    for beta in vcat(1e-3:1e-3:9e-3,1e-2:1e-2:1e-1)#keys(data[gname])
        #single graph single beta get statistic
        agg_data = data[gname][beta]

        push!(x_gname,qimpact_metric(agg_data))
        push!(y_gname,spectral_metric(pl_lam,pl_deg,raw_lam,raw_deg))
    end
    push!(xs,x_gname)
    push!(ys,y_gname)
end


xs = hcat(xs...)
ys = hcat(ys...)

#highlight beta=0.03 to see other plots 
f = scatter(xs,ys,leg=false,xscale=:log10,markerstrokewidth=0,
    xlims=(0.5,max(35,1.2*maximum(xs)+5)),ylims=(0.53,0.68))

scatter!(f,xs[12,1:end],ys[12,1:end],markerstrokewidth=0,color=:black,
    markershape=:star4, markersize=8)


# annotate!(f,150,0.68,text("5", 8))
# annotate!(f,150,0.65,text("10", 8))
# annotate!(f,90,0.643,text("15", 8))
# annotate!(f,35,0.625,text("20", 8))
# annotate!(f,40,0.615,text("25", 8))
# annotate!(f,40,0.568,text("30", 8))
# annotate!(f,28,0.558,text("35", 8))
# annotate!(f,23,0.55,text("40", 8))
# annotate!(f,13,0.57,text("45", 8))
# annotate!(f,18,0.583,text("50", 8))

annotate!(f,0.8,0.675,text("5", 8))
annotate!(f,0.8,0.655,text("10", 8))
annotate!(f,0.9,0.643,text("15", 8))
annotate!(f,0.8,0.625,text("20", 8))
annotate!(f,0.8,0.615,text("25", 8))
annotate!(f,0.8,0.566,text("30", 8))
annotate!(f,0.8,0.558,text("35", 8))
annotate!(f,0.8,0.55,text("40", 8))
annotate!(f,0.8,0.572,text("45", 8))
annotate!(f,0.8,0.583,text("50", 8))


plot!(f,[1;1],[ylims(f)...])
plot!(f,[1;xlims(f)[2]],[1;1],fillrange=[0; 0],
    alpha=0.15,color=3)

xlabel!(f,"Total Quarantine Infs Ratio (PL/Raw) ")
ylabel!(f,"Spectral Ratio (PL/Raw) ")

# ----------------------------------------------------------
# study plot
gnames = ["study-11-2022-1.smat","study-11-2022-10.smat",
    "study-11-2022-20.smat","study-11-2022-30.smat","study-11-2022-40.smat",
    "study-11-2022-50.smat"]

data = get_rewired_data(gnames)

xs = Vector{Vector}()
ys = Vector{Vector}()
lams = readdlm("pipeline/data/dominant-eigvals-0.txt")

for gname in gnames
    println("working on graph $gname")
    x_gname= Vector{Float64}()
    y_gname= Vector{Float64}()

    #load graph and get spectral information 
    A = loadGraph(gname,"input/graphs/")
    raw_deg = vec(sum(A;dims=2))
    raw_lam = lams[findfirst(lams[:,1].==gname),2]

    #PL rewired graph 
    A = loadGraph("rewired-10000.0-$gname","pipeline/graphs/")
    pl_deg = vec(sum(A;dims=2))
    pl_lam = lams[findfirst(lams[:,1].=="rewired-10000.0-$(gname)"),2]  
    
    for beta in vcat(1e-3:1e-3:9e-3,1e-2:1e-2:1e-1)#keys(data[gname])
        #single graph single beta get statistic
        agg_data = data[gname][beta]

        push!(x_gname,qimpact_metric(agg_data))
        push!(y_gname,spectral_metric(pl_lam,pl_deg,raw_lam,raw_deg))
    end
    push!(xs,x_gname)
    push!(ys,y_gname)
end


xs = hcat(xs...)
ys = hcat(ys...)

#highlight beta=0.03 to see other plots 
f = scatter(xs,ys,leg=false,xscale=:log10,markerstrokewidth=0,
    xlims=(0.5,max(35,1.2*maximum(xs)+5)),ylims=(0.53,0.8))

scatter!(f,xs[12,1:end],ys[12,1:end],markerstrokewidth=0,color=:black,
    markershape=:star4, markersize=8)

annotate!(f,0.8,0.79,text("initial", 8))
annotate!(f,0.8,0.655,text("10", 8))
annotate!(f,0.8,0.625,text("20", 8))
annotate!(f,0.8,0.566,text("30", 8))
annotate!(f,0.8,0.55,text("40", 8))
annotate!(f,0.8,0.583,text("50", 8))


plot!(f,[1;1],[ylims(f)...],color=3)
plot!(f,[1;xlims(f)[2]],[1;1],fillrange=[0; 0],
    alpha=0.15,color=3)

xlabel!(f,"Total Quarantine Infs Ratio (PL/Raw) ")
ylabel!(f,"Spectral Ratio (PL/Raw) ")

plot!(f,dpi=500)
Plots.savefig(f,"code/paper-figs/eigenvalue-plots/eigval-plot-spatial-plot.png")

