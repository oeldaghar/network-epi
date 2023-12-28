#after running diffusions, this script generates the aggreagted plotting data in the 
# future, can expand to all plotting data 

#store aggreagted data in /pipeline/data/[graphname]/diffusions/[diffusion type]/plotting-data
#load data in, aggregate and save to the above directories depending on graph 

using MatrixNetworks
using SparseArrays
using LinearAlgebra
using Distributions
using DelimitedFiles
using ProgressMeter
using Random

#allows for execution from command line as well as an ide so files can be run modularlly
mainDir = joinpath("/",split(abspath(""),"/")[1:findlast("network-epi" .== split(abspath(""),"/"))]...)

include(joinpath(mainDir,"code","graph-io.jl"))
include(joinpath(mainDir,"code","data-io.jl"))

ENV["GKSwstype"] = "100"

#---------------------------------------------------------------
#plotting functions 
using Plots
using Measures
function qpercent_heatmaps(gname::String;
        betas::Vector{Float64}=vcat(1e-3:1e-3:1e-2,2e-2:1e-2:1e-1),
        gammas::Vector{Float64}=[0.05],
        gpath::String="pipeline/graphs/",
        dloc::String="pipeline/data/",
        dtype::String="tinfs",
        method::String="seir",
        rewiring_type1::String="rewired",
        rewiring_type2::String="er",
        diffusion_type::String="uniform",
        ntrials::Int=50,
        qpercents::Vector{Int}=collect(0:15))


    pyplot()
    #for saving figs
    full,aggregated = Vector(),Vector()

    #main loop 
    @showprogress for beta in betas
        for gamma in gammas
            #should be a more julia centric way to do this.
            data,aggregated_data,clabel,rewiring_fractions1,rewiring_fractions2 = get_plotting_data(gname,beta=beta,
                gamma=gamma,gpath=gpath,dloc=dloc,dtype=dtype,method=method,
                rewiring_type1=rewiring_type1,rewiring_type2=rewiring_type2,
                diffusion_type=diffusion_type,ntrials=ntrials,qpercents=qpercents)
            #rewiring percentages
            rps = vcat(rewiring_fractions1,["0.0"],rewiring_fractions2)
            #full heatmap
            f = Plots.heatmap(data,
                ylabel="quarantine percentage",
                xlabel="     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent",
                title="$(gname[1:end-5]) - $(uppercase(method))($beta,$gamma)\n$(ntrials) node sample",
                colorbar_title=clabel,c=:heat,clims=(0,1),
                yticks=(collect(round(Int,ntrials/2):ntrials:size(data,1)),qpercents),xrotation = 90)

            Plots.plot!(f,xticks=(1:length(rps), vec(rps)),margins=9Measures.mm)
            push!(full,f)

            #heatmap using aggregated diffusion data.
            f = Plots.heatmap(aggregated_data,
                ylabel="quarantine percentage",
                xlabel="     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent",
                title="$(gname[1:end-5]) - $(uppercase(method))($beta,$gamma)\n$(ntrials) node sample",
                colorbar_title=clabel,c=:heat,clims=(0,1),
                yticks=(0.5:length(qpercents)-0.5,qpercents),xrotation = 90)

            Plots.plot!(f,xticks=(1:length(rps), vec(rps)),margins=9Measures.mm)
            push!(aggregated,f)
        end
    end    
    return full,aggregated
end


#gname,betas -> read all data from files and prep for plotting

#given single gname and beta, load diffusion data


function load_beta_data(gname::String;
          method::String="seir")

    g = canonical_graph_name(gname)
    #get params for loading files 
    betas = get_betas(g)
    sort!(betas)

    #line up gnames same as how we plot them 
    rps = sort_fnames(g,rewiring_type="rewired")[2]

    gnames = map(x->"rewired-$x-$gname",reverse(rps))
    gnames = vcat(gnames,g)
    gnames = vcat(gnames,map(x->"er-$x-$gname",rps))
    
    #load in data   
    data = Dict{Float64,Matrix}() #Dict{String,Dict}() #gname -> beta -> total
    dloc = "pipeline/data/$(g[1:end-5])/diffusions/uniform/"

    @showprogress for beta in betas
        tmp = Vector()
        for h in gnames
            tinfs = sum.(read_inf_data(h,dloc=dloc,beta=beta,dtype="tinfs",method=method))
            push!(tmp,tinfs)
        end
        data[beta] = hcat(tmp...)
    end

    return data,rps
end

function make_beta_rewiring_data(beta_data::Dict{Float64,Matrix},qpercent::Union{Float64,Int}=0)
    @assert(qpercent in collect(0:15),"qpercent outside valid range")
    offset = qpercent*50 #ntrials = 50
    betas = sort(collect(keys(beta_data)))
    result = vcat(map(x->beta_data[x][offset.+(1:50),:],betas)...)
    return result 
end


#TODO normalize by number of nodes in base graph but only load in graph once

# #for v1.8.0 plotting 
function qpercent_contours(gname::String;
        betas::Vector{Float64}=vcat(1e-3:1e-3:1e-2,2e-2:1e-2:1e-1),
        gammas::Vector{Float64}=[0.05],
        gpath::String="pipeline/graphs/",
        dloc::String="pipeline/data/",
        dtype::String="tinfs",
        method::String="seir",
        rewiring_type1::String="rewired",
        rewiring_type2::String="er",
        diffusion_type::String="uniform",
        ntrials::Int=50,
        qpercents::Vector{Int}=collect(0:15),
        add_labels::Bool=true,
        colorbar::Bool=true)


    # pyplot()
    #for saving figs
    full,aggregated = Vector(), Vector()

    #main loop 
    @showprogress for beta in betas
        for gamma in gammas
            #should be a more julia centric way to do this.
            data,aggregated_data,clabel,rewiring_fractions1,rewiring_fractions2 = get_plotting_data(gname,beta=beta,
                gamma=gamma,gpath=gpath,dloc=dloc,dtype=dtype,method=method,
                rewiring_type1=rewiring_type1,rewiring_type2=rewiring_type2,
                diffusion_type=diffusion_type,ntrials=ntrials,qpercents=qpercents)
            #rewiring percentages
            rps = vcat(rewiring_fractions1,["0.0"],rewiring_fractions2)
            # #contour plot for full data
            # xlab,ylab,title = "","",""
            # if add_labels
            #     ylab = "quarantine percentage"
            #     xlab = "     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent"
            #     title = "$(gname[1:end-5]) - $(uppercase(method))($beta,$gamma)\n$(ntrials) node sample"
            # end

            f = Plots.contourf(1:length(rps),1:size(data,1),(x,y)->data[y,x],
                # ylabel=ylab,
                # xlabel=xlab,
                # title=title,
                # colorbar_title=clabel,
                c=:heat,clims=(0,1),
                # yticks=(collect(round(Int,ntrials/2):ntrials:size(data,1)),qpercents),
                # xrotation = 90,
                levels=0.0:0.2:1.0)
                
            # Plots.plot!(f,xticks=(1:length(rps), vec(rps)),margins=9Measures.mm)
            push!(full,f)
         
            #contour plot for aggreagted data
            f1 = Plots.contourf(1:length(rps),1:length(qpercents),(x,y)->aggregated_data[y,x],
                # ylabel=ylab,
                # xlabel=xlab,
                # title=title,
                # colorbar=colorbar,
                c=:heat,
                # colorbar_title=clabel,
                clims=(0,1),
                # yticks=(1:length(qpercents),qpercents),
                # xrotation = 90,
                levels=0.0:0.2:1.0)
                
            # Plots.plot!(f1,xticks=(1:length(rps), vec(rps)),margins=9Measures.mm)
            push!(aggregated,f1)

        end
    end    
    return full,aggregated
end

function qpercent_contours_test(gname::String;
    betas::Vector{Float64}=vcat(1e-3:1e-3:1e-2,2e-2:1e-2:1e-1),
    gammas::Vector{Float64}=[0.05],
    gpath::String="pipeline/graphs/",
    dloc::String="pipeline/data/",
    dtype::String="tinfs",
    method::String="seir",
    rewiring_type1::String="rewired",
    rewiring_type2::String="er",
    diffusion_type::String="uniform",
    ntrials::Int=50,
    qpercents::Vector{Int}=collect(0:15),
    add_labels::Bool=true,
    colorbar::Bool=true)


    # pyplot()
    #for saving figs
    full,aggregated = Vector(), Vector()

    #main loop 
    @showprogress for beta in betas
        for gamma in gammas
            #should be a more julia centric way to do this.
            data,aggregated_data,clabel,rewiring_fractions1,rewiring_fractions2 = get_plotting_data(gname,beta=beta,
                gamma=gamma,gpath=gpath,dloc=dloc,dtype=dtype,method=method,
                rewiring_type1=rewiring_type1,rewiring_type2=rewiring_type2,
                diffusion_type=diffusion_type,ntrials=ntrials,qpercents=qpercents)
            #rewiring percentages
            rps = vcat(rewiring_fractions1,["0.0"],rewiring_fractions2)
            # #contour plot for full data
            # xlab,ylab,title = "","",""
            # if add_labels
            #     ylab = "quarantine percentage"
            #     xlab = "     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent"
            #     title = "$(gname[1:end-5]) - $(uppercase(method))($beta,$gamma)\n$(ntrials) node sample"
            # end

            f = Plots.contourf(1:length(rps),1:size(data,1),(x,y)->data[y,x],
                # ylabel=ylab,
                # xlabel=xlab,
                # title=title,
                # colorbar_title=clabel,
                c=:heat,clims=(0,1),
                # yticks=(collect(round(Int,ntrials/2):ntrials:size(data,1)),qpercents),
                # xrotation = 90,
                levels=0.0:0.05:1.0)
                
            # Plots.plot!(f,xticks=(1:length(rps), vec(rps)),margins=9Measures.mm)
            push!(full,f)
        
            #contour plot for aggreagted data
            f1 = Plots.contourf(1:length(rps),1:length(qpercents),(x,y)->aggregated_data[y,x],
                # ylabel=ylab,
                # xlabel=xlab,
                # title=title,
                # colorbar=colorbar,
                c=:heat,
                # colorbar_title=clabel,
                clims=(0,1),
                # yticks=(1:length(qpercents),qpercents),
                # xrotation = 90,
                levels=0.0:0.05:1.0)
                
            # Plots.plot!(f1,xticks=(1:length(rps), vec(rps)),margins=9Measures.mm)
            push!(aggregated,f1)

        end
    end    
    return full,aggregated
end

#for plotting triangle diffusion plots 
function qpercent_contours_triangles_diffusions(gname::String;
        betas::Vector{Float64}=vcat(1e-3:1e-3:1e-2,2e-2:1e-2:1e-1),
        gammas::Vector{Float64}=[0.05],
        gpath::String="pipeline/graphs/",
        dloc::String="pipeline/data/",
        dtype::String="tinfs",
        method::String="seir",
        ntrials::Int=50,
        qpercents::Vector{Int}=collect(0:15),
        add_labels::Bool=true,
        colorbar::Bool=true)


    #for saving figs
    full,aggregated = Vector(), Vector()

    #main loop 
    @showprogress for beta in betas
        for gamma in gammas
            #should be a more julia centric way to do this.
            data,aggregated_data,clabel,rewiring_fractions = get_triangle_diffusions_plotting_data(gname,beta,gamma,
                                        dloc,dtype,method,ntrials)
            
            #rewiring percentages
            rps = rewiring_fractions

            f = Plots.contourf(1:length(rps),1:size(data,1),(x,y)->data[y,x],
                # ylabel=ylab,
                # xlabel=xlab,
                title="$gname\n$method($beta,$gamma)",
                colorbar_title=clabel,
                c=:heat,clims=(0,1),
                yticks=(collect(round(Int,ntrials/2):ntrials:size(data,1)),qpercents),
                # xrotation = 90,
                levels=0.0:0.2:1.0)
                
            # Plots.plot!(f,xticks=(1:length(rps), vec(rps)),margins=9Measures.mm)
            push!(full,f)
         
            #contour plot for aggreagted data
            f1 = Plots.contourf(1:length(rps),1:length(qpercents),(x,y)->aggregated_data[y,x],
                # ylabel=ylab,
                # xlabel=xlab,
                # title="$gname\n$method($beta,$gamma)",
                # colorbar_title=clabel,
                c=:heat,clims=(0,1),
                # yticks=(1:length(qpercents),qpercents),
                # xrotation = 90,
                levels=0.0:0.2:1.0)
                
            # Plots.plot!(f1,xticks=(1:length(rps), vec(rps)),margins=5Measures.mm)
            push!(aggregated,f1)

        end
    end    
    return full,aggregated
end


#eigenvalues
function _align_rewiring_percents(rps::Vector{Float64};
            full_rps::Vector{Float64}=vec([0.01 0.05 0.1 0.5 1.0 5.0 10.0 25.0 50.0 100.0 1000.0 10000.0]))

    full_rps = vcat(reverse(full_rps), [0], full_rps)    
    @assert(all(x in full_rps for x in rps))

    result = Vector{Int}()

    i,j = 1,1
    while (i<=lastindex(rps) && j<=lastindex(full_rps))
        if rps[i]!=full_rps[j]
            j+=1
        else
            push!(result,j)
            i+=1
        end
    end
    return result
end

# _align_rewiring_percents([10.0; 1.0; 0; 1.0; 10.0])
# [6; 8; 13; 18; 20]


function plot_eigenvalues(gname::String;gpath::String="pipeline/graphs/",
                full_rps::Vector{Float64}=vec([0.01 0.05 0.1 0.5 1.0 5.0 10.0 25.0 50.0 100.0 1000.0 10000.0]),
                test::Bool=false)
    
    full_rps = vcat(reverse(full_rps), [0], full_rps)    

    # pyplot()

    println("working on graph $gname")
    gnames = getgnames(gname,gpath)
    gname = gnames[1]
    filter!(c->startswith(c,"er-") || startswith(c,"rewired-"),gnames)
    
    if !startswith(gname,"cn-") #toss out cases like Penn picking up cn-Penn
        filter!(x->!occursin("cn-",x),gnames)
    end

    @assert(length(gnames)>0)
    rewiring_ps = sort(unique(parse.(Float64,map(x->split(x,"-")[2],gnames))))

    #load eigvalue data 
    eig_data = readdlm("pipeline/data/dominant-eigvals-0.txt")

    graph_eigs = Vector{Float64}()
    
    for ps in reverse(rewiring_ps)
        tmp_gname = "rewired-$ps-$gname"
        push!(graph_eigs,eig_data[findfirst(eig_data[:,1].==tmp_gname),2])
    end
    
    push!(graph_eigs,eig_data[findfirst(eig_data[:,1].==gname),2])

    for ps in rewiring_ps
        tmp_gname = "er-$ps-$gname"
        push!(graph_eigs,eig_data[findfirst(eig_data[:,1].==tmp_gname),2])
    end
    
    #plot the data 
    rps = vcat(reverse(rewiring_ps), [0], rewiring_ps)
    rps_inds = collect(eachindex(full_rps))
    if rps!=full_rps #one of the larger graphs, didn't do full rewiring range
        #align indices of reduced data to indices of full data
        rps_inds = _align_rewiring_percents(rps)
    end 
    f = Plots.scatter(rps_inds,graph_eigs,
        # ylabel="",
        # ylabel="Dominant Eigenvale",
        # xlabel="     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent",
        # xlabel="",
        # title="$(gname[1:end-5])",
        xrotation = 90, markerstrokewidth=0,markersize=5)
    
    # return f 
    ind = Int(findfirst(rps.==0))
    Plots.scatter!(f,[rps_inds[ind]],[graph_eigs[ind]],markerstrokewidth=0,
            markersize=9)
    if test #add labels + xticks
        Plots.plot!(f,xticks=(rps_inds, string.(vec(rps))),margins=9Measures.mm)
                        # ylabel="Dominant Eigenvale",title="$(gname[1:end-5])")
    end
    Plots.plot!(f,leg=false)

    return f
end


function plot_eigenvalues(gnames::Vector{String};gpath::String="pipeline/graphs/")
    
    pyplot()
    gnames = getgnames.(gnames,gpath)
    #load eigvalue data 
    eig_data = readdlm("pipeline/data/dominant-eigvals-0.txt")

    data = Vector{Vector{Float64}}()
    rewiring_percents = Vector{Vector{Float64}}()
    
    for i in eachindex(gnames)
        tmp_gnames = gnames[i]
        
        tmp_gname = tmp_gnames[1]
        @show tmp_gname
        filter!(c->startswith(c,"er-") || startswith(c,"rewired-"),tmp_gnames)
        
        if !startswith(tmp_gname,"cn-") #toss out cases like Penn picking up cn-Penn
            filter!(x->!occursin("cn-",x),tmp_gnames)
        end

        @assert(length(tmp_gnames)>0)
        rewiring_ps = sort(unique(parse.(Float64,map(x->split(x,"-")[2],tmp_gnames))))

        
        graph_eigs = Vector{Float64}()
        @show findfirst(eig_data[:,1].==tmp_gname)
        for ps in reverse(rewiring_ps)
            gname = "rewired-$ps-$tmp_gname"
            push!(graph_eigs,eig_data[findfirst(eig_data[:,1].==gname),2])
        end
        
        push!(graph_eigs,eig_data[findfirst(eig_data[:,1].==tmp_gname),2])

        for ps in rewiring_ps
            gname = "er-$ps-$tmp_gname"
            push!(graph_eigs,eig_data[findfirst(eig_data[:,1].==gname),2])
        end
        push!(data,graph_eigs)
        push!(rewiring_percents,vcat(reverse(rewiring_ps), [0], rewiring_ps))
    end
    

    #plot the data 
    return data
    # rps = vcat(reverse(rewiring_ps), [0], rewiring_ps)
    # f = Plots.scatter(1:length(rps),graph_eigs,
    #     ylabel="Dominant Eigenvale",
    #     xlabel="     Powerlaw ⟺ Original ⟺ Erdos Renyi\nRewiring Percent",
    #     title="$(gname[1:end-5])",
    #     xrotation = 90, markerstrokewidth=0)
    
    # # return f 
    # ind = Int(findfirst(rps.==0))
    # Plots.scatter!(f,[ind],[graph_eigs[ind]],markerstrokewidth=0,
    #         markersize=10)
    # Plots.plot!(f,xticks=(1:length(rps), string.(vec(rps))),margins=9Measures.mm, leg = false)
    #
    # return f
end
 

#graphs to use and their position
hs = ["commutes-all","mexico", "filtered", #row1
    "cn-moduillinois", "cn-Penn", "cn-modWiscon", #row2...
    "dblp","enron","anon", #row 3
    "cit-HepPh", "slashdot", "flickr", 
    "geometric","study-25-150","cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5.smat"]#"cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.9-0.0-100-5.smat"]
    #  "rewired-10000.0-modmexico","rewired-10000.0-cn-moduillinois","er-10000.0-dblp-cc"]
gnames = map(h->h=="" ? "" : getgnames(h,"pipeline/graphs/")[1], hs) 
titles = ["US Commutes", "Mexico City\nTrace", "Filtered\nUS Flows",
    "Sparsified\nCollege-Illinois", "Sparsified\nCollege-Penn", "Sparsified\nCollege-Wisc.",    
    "Collaboration", "Email", "Facebook\nInteractions",
    "Citation", "Slashdot","Flickr",    
    "Local\nGeometric", "Geometric\nCommunities", "Random Walk\nCommunities"] 
ybounds = [(91,280),(1,500),(18,65),
    (3,31),(2.5,38),(3,25),
    (2,82),(4,128),(4,23),
    (22,82),(7,142),(-45,1350),    
    (14.8,18.6),(29,105),(5,31),
]

ps = map(x -> plot_eigenvalues(x),gnames)
ps = vcat(ps...)

figs = deepcopy(ps)
#remove text, and retitle
for (ind,title) in enumerate(titles)
    f = figs[ind]
    # plot!(f,grid=false,xlabel="",ylabel="",title=title,
    #     titlefontsize=12,margins=0Measures.mm,xaxis=false)
    plot!(f,grid=false,xlabel="",ylabel="",xaxis=false,
        leftmargin=1Measures.mm,
        bottommargin=-5Measures.mm)

    ymin,ymax = ylims(f)
    annotate!(f, (mean(xlims(f)), ymin+0.25(ymax-ymin), text(title, 12, :left, RGB(0.0,0.0,0.0), :center)))

    #increase limit for y to avoid cutting off markers
    plot!(f,ylims = ybounds[ind])
end


#labels for specific plots
using LaTeXStrings
plot!(figs[13],ylabel=L"\lambda_1(A)",leftmargin=5Measures.mm,yguidefontsize=20)
# plot!(figs[14],xlabel="Rewiring Fraction",bottommargin=0Measures.mm,xguidefontsize=20,
#     xaxis=true,)

#plotting
# newf = Plots.plot(figs...,
#         layout=grid(5,3,heights=[0.2, 0.2, 0.2, 0.2, 0.2],widths=[0.33,0.33,0.33]/(0.33+0.33+0.33)),
#         size=(1000,1200),
#         link=:x)
heights = [1, 1, 1, 1, 1.0]
heights./=sum(heights)
widths = [1.0,1,1]
widths./=sum(widths)

newf = Plots.plot(figs...,
    layout=grid(5,3,heights=heights,widths=widths),
    size=(1000,600),
    link=:x,
    top_margin=0Measures.mm,
    bottom_margin=-5Measures.mm,
    right_margin=0Measures.mm)


# plot!(newf,inset = (bbox(0, 0, 0.05, 0.05, :center,:center)))
# annotate!(newf,inset = (bbox(0, 0, 0.05, 0.05, :center,:center)))
# annotate!(newf, (mean(xlims(f)), ymin+0.25(ymax-ymin), text(title, 12, :left, RGB(0.0,0.0,0.0), :center)))
plot!(newf,dpi=1200)
Plots.savefig(newf,"code/paper-figs/eigenvalue-plots/eigvals-plot-base.png")
Plots.savefig(newf,"code/paper-figs/eigenvalue-plots/eigvals-plot-base.pdf")


inds = [1;2;3;13;14;4;7;9;10;12]
newf = Plots.plot(figs[inds]...,
    layout=grid(2,5,heights=[0.5,0.5],widths=[1/5 for i=1:5]),
    size=(1500,350),
    link=:x,
    top_margin=1Measures.mm,
    bottom_margin=-1Measures.mm,
    right_margin=0Measures.mm,
    ygrid=false)
    # grid=(:y, :black, :dot, 1, 0.9))
plot!(newf,dpi=800)
Plots.savefig(newf,"code/paper-figs/eigenvalue-plots/eigvals-plot-v2.png")
    


gname = getgnames("draft","input/graphs/")[1]
ps = qpercent_contours(gname,betas=get_betas(gname),colorbar=false,add_labels=false)[2]
ps = vcat(ps...)

f = deepcopy(ps[10])
plot!(f,yticks=(1.0:16.0,0:15))
rps = ["10000";"1000";"100";"50";"25";"10";"5";"1";"0.5";"0.1";"0.05";"0.01";"0";"0.01";"0.05";"0.1";"0.5";"1";"5";"10";"25";"50";"100";"1000";"10000"]
plot!(f,xticks=(1:25,rps),xrotation=90,guidefontweight="bold")
plot!(f,colorbar=false)
plot!(f,dpi=1000)
Plots.savefig(f,"code/paper-figs/heatmaps/$(gname[1:end-5])-uplot.png")


#########################################################################################################
####################### UPLOTS #####################################
####################################################################

#after testing we've found the betas we want to use 
#graphs to use and their position
hs = ["commutes-all","mexico", "filtered", #row1
    "cn-moduillinois", "cn-Penn", "cn-modWiscon", #row2...
    "dblp","enron","anon", #row 3
    "cit-HepPh", "slashdot", "flickr", 
    "geometric","study-25-150",#"cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.9-0.0-100-5.smat",
    "cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5.smat"]
    #  "rewired-10000.0-modmexico","rewired-10000.0-cn-moduillinois","er-10000.0-dblp-cc"]
gnames = map(h->h=="" ? "" : getgnames(h,"pipeline/graphs/")[1], hs) 
titles = ["US\nCommutes", "Mexico City\nTrace", "Filtered\nUS Flows",
    "Sparsified\nCollege-Illinois", "Sparsified\nCollege-Penn", "Sparsified\nCollege-Wisc.",    
    "Collaboration", "Email", "Facebook\nInteractions",
    "Citation", "Slashdot","Flickr",    
    "Local\nGeometric", "Geometric\nCommunities", "Random Walk\nCommunities"] 

# betas = [[0.003],[0.05],[0.03],
#     [0.1],[0.1],[0.1],
#     [0.1],[0.06],[0.1],
#     [0.008],[0.04],[0.04],
#     [0.06],[0.03],[0.15]]

betas = [[0.003],[0.05],[0.03],
    [0.1],[0.1],[0.1],
    [0.1],[0.06],[0.1],
    [0.008],[0.04],[0.04],
    [0.06],[0.02],[0.175]]

ps = map(x -> qpercent_contours(x[1],betas=x[2],colorbar=false,add_labels=false)[2],zip(gnames,betas))
ps = vcat(ps...)


#for tinkering and making adjustments without having to remake figs from scratch
figs = deepcopy(ps)

#remove text, and retitle
#this is specific to Julia v1.5.1 w/ pyplot backend. not sure why there's issues..
for (ind,title) in enumerate(titles)
    f = figs[ind]
    # title = title*"\nSEIR($(betas[ind][1]),0.05)"
    # title = title*"\nSEIR($(betas[ind][1]),0.05)"
    # plot!(f,framestyle=:none,xlabel="",ylabel="",title=title,
    #     titlefontsize=12,margins=0Measures.mm)
    ymin,ymax = ylims(f)
    plot!(f,framestyle=:none,xlabel="",ylabel="",title="",
        xguidefontsize=12,
        right_margin=0Measures.mm,
        colorbar=false)
    annotate!(f,[(mean(xlims(f)),ymin+0.8*(ymax-ymin),text(title, 16, :center, :center,RGB(0.0,0.0,0.0)))])
end

function add_tickmarks(f,xmirror=false,ymirror=false)
    plot!(f,xmirror=xmirror,ymirror=ymirror,
        tick_direction=:out,
        xticks = (range(1.5,24.5,25),["" for i=1:25]),xrotation=90,
        top_margin=2Measures.mm,
        yticks = (range(1.5,15.5,16),["" for i=1:16]),
        colorbar = false,framestyle=true,
        xlims=(1.05,24.95),ylims=(1.05,15.95))
    return f
end

#modifying exterior plots 
f = deepcopy(figs[1])
for (i,f) in enumerate(figs)
    xmirror,ymirror = false,false
    if i in [1;2;3]
        xmirror = true
    end
    if i in [3;6;9;12;15]
        ymirror = true
    end
    f = add_tickmarks(f,xmirror,ymirror)
    
    if i in 4:11
        plot!(f,xticks=false)
    end

    if i in 2:3:15
        plot!(f,yticks=false)
    end
    figs[i] = f
end

#for saving individual plots
# for f in figs
#     plot!(f,xticks=false,yticks=false)
# end
# figs[1]

#manually handle flickr
f = figs[12]
# plot!(f,xticks=false,xlims=(1,11))
plot!(f,xticks=(range(1.5,10.5,11),["" for i=1:11]), tick_direction=:out,xlims=(1,11))
f


#plotting
newf = Plots.plot(figs...,
        layout=grid(5,3,heights=[0.2, 0.2, 0.2, 0.2, 0.2],widths=[0.33,0.33,0.33]/(0.33+0.33+0.33)),
        size=(1000,1200),
        # margins=0Measures.mm,
        bottom_margin=-2*Measures.mm,
        right_margin=-1*Measures.mm,
        left_margin=-2*Measures.mm,
        top_margin=-0.5*Measures.mm,
        link=:y)
plot!(newf,dpi=1500)
Plots.savefig(newf,"code/paper-figs/heatmaps/uplot-final.png")
Plots.savefig(newf,"code/paper-figs/heatmaps/uplot-final.pdf")

#high resolution image of mexico city for explatory figure and for colorbar
f = deepcopy(ps[2])
plot!(f,yticks=(1.0:16.0,0:15))
rps = ["10000";"1000";"100";"50";"25";"10";"5";"1";"0.5";"0.1";"0.05";"0.01";"0";"0.01";"0.05";"0.1";"0.5";"1";"5";"10";"25";"50";"100";"1000";"10000"]
plot!(f,xticks=(1:25,rps),xrotation=90,guidefontweight="bold")
plot!(f,colorbar=false)
plot!(f,dpi=1000)
Plots.savefig(f,"code/paper-figs/heatmaps/mexico-city-uplot.png")

#high resolution image of mexico city for explatory figure and for colorbar
f = deepcopy(figs[2])
plot!(f,
        top_margin = 1Measures.mm, 
        bottom_margin = -1*Measures.mm,
        right_margin = 2*Measures.mm, 
        left_margin = -9*Measures.mm,
        colorbar=true)
plot!(f,dpi=2000)
Plots.savefig(f,"code/paper-figs/heatmaps/uplot-colorbar.png")


#for saving individual plots
for (ind,f) in enumerate(figs)
    f = add_tickmarks(f)
    plot!(f,dpi=800)
    savefig(f,"code/paper-figs/heatmaps/uplot-standalone-$(gnames[ind][1:end-5]).png")
end



### alternate layout 
figs = deepcopy(ps)
for (ind,title) in enumerate(titles)
    f = figs[ind]
    # title = title*"\nSEIR($(betas[ind][1]),0.05)"
    # title = title*"\nSEIR($(betas[ind][1]),0.05)"
    # plot!(f,framestyle=:none,xlabel="",ylabel="",title=title,
    #     titlefontsize=12,margins=0Measures.mm)
    ymin,ymax = ylims(f)
    plot!(f,framestyle=:none,xlabel="",ylabel="",title="",
        xguidefontsize=12,
        # bottommargin=Measures.mm,
        bottom_margin=-8*Measures.mm,
        right_margin=-5*Measures.mm,
        left_margin=-10*Measures.mm,
        top_margin=-5*Measures.mm,
        colorbar=false)
    annotate!(f,[(mean(xlims(f)),ymin+0.7*(ymax-ymin),text(title, 16, :center, :center,RGB(0.0,0.0,0.0)))])
end
figs[end]


###alternate layout
inds = [1;2;3;13;14;4;7;9;10;12]
newfigs = deepcopy(figs[inds])
for (i,f) in enumerate(newfigs)
    xmirror,ymirror = false,false
    if i in [1;2;3;4;5]
        xmirror = true
    end
    if i in [5;10]
        ymirror = true
    end
    f = add_tickmarks(f,xmirror,ymirror)
    
    if i in vcat(2:4,7:9)
        plot!(f,yticks=false)
    end
    newfigs[i] = f
end

#manually handle flickr
f = newfigs[10]
# plot!(f,xticks=false,xlims=(1,11))
plot!(f,xticks=(range(1.5,10.5,11),["" for i=1:11]), tick_direction=:out,xlims=(1,11))
f

newf = Plots.plot(newfigs...,
        layout=grid(2,5,heights=[0.5,0.5],widths=[1/5 for i=1:5]),
        size=(1500,500),
        bottom_margin=-3Measures.mm,
        left_margin=-2Measures.mm,
        right_margin=-1Measures.mm,
        link=:y,
        dpi=500)

Plots.savefig(newf,"code/paper-figs/heatmaps/uplot-v2.png")



#### GENERATING PLOTS FOR LFR + SPATIAL NETWORKS FIGURES

#LFR 
gnames = getgnames("cl-lfr","input/graphs/")[1:end-1]
figs = []
for gname in gnames
    fs = qpercent_contours(gname,betas=[0.175])
    push!(figs,fs[2]...)
end
for (i,f) in enumerate(figs)
    gname = gnames[i]
    plot!(f,colorbar=false,framestyle=:none,
        margins=-10Measures.mm,dpi=500)
    Plots.savefig(f,"code/paper-figs/heatmaps/uplot-$gname.png")
end


#Spatial 
# gnames = ["study-11-2023-0-noweak.smat","study-11-2022-1.smat","study-11-2022-10.smat",
#     "study-11-2022-20.smat","study-11-2022-30.smat",
#     "study-11-2022-40.smat","study-11-2022-50.smat"]

gnames = [
    "study-11-2023-0-noweak.smat",
    "study-11-2022-1.smat",
    "study-11-2023-1-longrange-1.smat",
    "study-11-2023-1-longrange-2.smat",
    "study-11-2023-1-longrange-3.smat",
    "study-11-2023-1-longrange-5.smat",
    "study-11-2023-1-longrange-8.smat",
    "study-11-2023-1-longrange-10.smat",
    ]
figs = []
for gname in gnames
    fs = qpercent_contours(gname,betas=[0.03])
    push!(figs,fs[2]...)
end

figs[1]
figs[2]
figs[3]
figs[4]
figs[5]
figs[6]
figs[7]

gname = "study-11-2022-1.smat"
gname = "study-11-2023-0-noweak.smat"
gname = "study-11-2023-1-longrange-15.smat"
fs = qpercent_contours(gname,betas=get_betas(gname))

ind = 15
get_betas(gname)[ind]
fs[2][ind]

findfirst(eig_data[:,1].==gname)

plot!.(figs,colorbar=false,framestyle=:none,
    left_margin=-7Measures.mm,
    right_margin=-1Measures.mm,
    bottom_margin=-5Measures.mm,
    top_margin=0Measures.mm,)


plot(figs...,layout=(1,5),size=(3200,600))
figs[1]
for (i,f) in enumerate(figs)
    gname = gnames[i]
    plot!(f,colorbar=false,framestyle=:none,
        margins=-10Measures.mm,dpi=500)
    Plots.savefig(f,"code/paper-figs/heatmaps/uplot-$(gname[1:end-5]).png")
end

#different parameters 
figs = []
tmpbeta = 5e-2
for gname in gnames
    fs = qpercent_contours(gname,betas=[tmpbeta])
    push!(figs,fs[2]...)
end
for (i,f) in enumerate(figs)
    gname = gnames[i]
    plot!(f,colorbar=false,framestyle=:none,
        margins=-10Measures.mm,dpi=500)
    Plots.savefig(f,"code/paper-figs/spatial-figs/uplot-$(gname[1:end-5])-seir-$tmpbeta.png")
end


# gnames[1]
# gnames[2]
# A = loadGraph(gnames[1],"input/graphs/")
# nnz(A)/size(A,1)
# B = loadGraph(gnames[2],"input/graphs/")
# nnz(B)/size(B,1)

# (nnz(A)/size(A,1))/(nnz(B)/size(B,1))*1e-1


##############################################################################
###################### TRIANGLE WEIGHTED UPLOTS ###########################
###########################################################################

hs = ["commutes-all","mexico", "filtered", #row1
    "cn-moduillinois", "cn-Penn", "cn-modWiscon", #row2...
    "dblp","enron","anon", #row 3
    "cit-HepPh", "slashdot", "flickr", 
    "geometric","study-25-150","cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5.smat"]
gnames = map(h->h=="" ? "" : getgnames(h,"pipeline/graphs/")[1], hs) 
titles = ["US\nCommutes", "Mexico City\nTrace", "Filtered\nUS Flows",
    "Sparsified\nCollege-Illinois", "Sparsified\nCollege-Penn", "Sparsified\nCollege-Wisc.",    
    "Collaboration", "Email", "Facebook\nInteractions",
    "Citation", "Slashdot","Flickr",    
    "Local\nGeometric", "Geometric\nCommunities", "Random Walk\nCommunities"] 

betas = [
    [0.0001],[0.003],[0.001],
    [0.18],[0.18],[0.18],
    [0.01],[0.001],[0.15],
    [0.001],[0.006],[0.03],
    [0.005],[5e-4],[0.06]]

#triangle weighted diffusion plots

ps = map(x -> qpercent_contours_triangles_diffusions(x[1],betas=x[2],colorbar=false,add_labels=false)[2],zip(gnames,betas))
ps = vcat(ps...)

annotation_locations = [
  (3.5,12), (3.5,12), (3,12) ,#row1
  (6,12), (6,12), (6,12), #row2
  (5,12), (4,12), (5,12), #row3
  (3,13), (3,13), (1.5,13),
  (4,12), (4,12), (8,12)
]

#for tinkering and making adjustments without having to remake figs from scratch
figs = deepcopy(ps)
#remove text, and retitle
#this is specific to Julia v1.5.1 w/ pyplot backend. not sure why there's issues..
for (ind,title) in enumerate(titles)
    f = figs[ind]
    # title = title*"\nSEIR($(betas[ind][1]),0.05)"
    # title = title*"\nSEIR($(betas[ind][1]),0.05)"
    # plot!(f,framestyle=:none,xlabel="",ylabel="",title=title,
    #     titlefontsize=12,margins=0Measures.mm)
    ymin,ymax = ylims(f)
    plot!(f,framestyle=:none,xlabel="",ylabel="",title="",
        xguidefontsize=12,
        bottom_margin=0*Measures.mm,
        top_margin = 0*Measures.mm,
        right_margin=0*Measures.mm,
        left_margin=0*Measures.mm,
        colorbar=false)
    xann,yann = annotation_locations[ind]
    annotate!(f,[(xann,yann,text(title, 16, :center, :center,RGB(0.0,0.0,0.0)))])
end

#adding tick marks 
function add_tickmarks_triangle(f,xmirror=false,ymirror=false)
    plot!(f,xmirror=xmirror,ymirror=ymirror,tick_direction=:out,
        xticks = (range(1.5,10.5,11),["" for i=1:11]),xrotation=90,
        top_margin=2Measures.mm,
        yticks = (range(1.5,15.5,16),["" for i=1:16]),
        tickfontsize=10,colorbar = false,framestyle=true,
        ylims=(1.05,15.95))
    return f
end

# figs[1]

#modifying exterior plots 
for (i,f) in enumerate(figs)
    xmirror,ymirror = false,false
    if i in [1;2;3]
        xmirror = true
    end
    if i in [3;6;9;12;15]
        ymirror = true
    end
    f = add_tickmarks_triangle(f,xmirror,ymirror)
    
    if i in 4:11
        plot!(f,xticks=false)
    end

    if i in 2:3:15
        plot!(f,yticks=false)
    end
    figs[i] = f
end
# figs[1]
#manually handle flickr
f = figs[12]
# plot!(f,xticks=false,xlims=(1,11))
plot!(f,xticks=(range(1.5,3.5,4),["" for i=1:4]), tick_direction=:out,xlims=(1,4))


#plotting
newf = Plots.plot(figs...,
        layout=grid(5,3,heights=[0.2, 0.2, 0.2, 0.2, 0.2],widths=[0.33,0.33,0.33]/(0.33+0.33+0.33)),
        size=(1000,1200),
        # margins=0Measures.mm,
        bottom_margin=-2*Measures.mm,
        right_margin=0*Measures.mm,
        left_margin=-2*Measures.mm,
        top_margin=0*Measures.mm,
        link=:y
        )
plot!(newf,dpi=1500)
Plots.savefig(newf,"code/paper-figs/heatmaps/uplot-triangles-v1.png")



##############################################
####### LOOKING AT NUMBER OF TRIANGLES #######
##############################################



hs = ["commutes-all","mexico", "filtered", #row1
    "cn-moduillinois", "cn-Penn", "cn-modWiscon", #row2...
    "dblp","enron","anon", #row 3
    "cit-HepPh", "slashdot", "flickr", 
    "geometric","study-25-150","cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5.smat"]
gnames = map(h->h=="" ? "" : getgnames(h,"pipeline/graphs/")[1], hs) 
titles = ["US\nCommutes", "Mexico City\nTrace", "Filtered\nUS Flows",
    "Sparsified\nCollege-Illinois", "Sparsified\nCollege-Penn", "Sparsified\nCollege-Wisc.",    
    "Collaboration", "Email", "Facebook\nInteractions",
    "Citation", "Slashdot","Flickr",    
    "Local\nGeometric", "Geometric\nCommunities", "Random Walk\nCommunities"] 

betas = [
    [0.0001],[0.003],[0.001],
    [0.18],[0.18],[0.18],
    [0.01],[0.001],[0.15],
    [0.001],[0.006],[0.03],
    [0.005],[5e-4],[0.06]]



#rewiring param versus triangle ratio to base graph 




##########################
# playing around with triangle rewiring 

using SparseArrays


"""
    commonNeighbors(A::SparseMatrixCSC)

returns a sparse matrix P with same structure as A with P[u,v] = number of triangles both u and v particpate in
"""
function commonNeighbors(A::SparseMatrixCSC)
    @assert issymmetric(A)
    P = deepcopy(A)
    fill!(nonzeros(P),0)
    rowval = rowvals(A)
    neighs = zeros(Bool, size(A,1))
    for j=1:size(A,1)
        # index neighbors
        for nzj in nzrange(A, j)
            neighs[rowval[nzj]] = true
        end
        # for each edge out of this node...
        for nzj in nzrange(A, j)
            i = rowval[nzj]
            score = 0
            for nzi in nzrange(A, i)
                w = rowval[nzi]
                if neighs[w] == true
                    # for each neighbor of i (w) that is
                    # also a neighbor of j neighs[w] = true
                    # increase the score
                    score += 1
                end
            end
                nonzeros(P)[nzj] = score
        end
        # reset indicies
        for nzj in nzrange(A, j)
            neighs[rowval[nzj]] = false
        end
    end
    return P
end

"""
    commonNeighbors(adjlist,u,v)

returns the number of neighbors shared by nodes u and v 
"""
function commonNeighbors(adjlist,u,v)
    neighs = zeros(Bool, lastindex(adjlist))
    score = 0
    for (node,_) in adjlist[u]
        neighs[node] = true
    end

    for (node,_) in adjlist[v]
        score+=neighs[node]
    end
    return score
end


"""
    sparse_to_adlist(A::SparseMatrixCSC)

converts a sparse matrix to adjacency list format where adjlist[u] = [(v,w_v)]
is the weighted adjacency list
"""
function sparse_to_adjlist(A::SparseMatrixCSC)::Vector{Vector{Tuple{Int,Int}}}
    @assert(issymmetric(A),"input matrix not symmetric")
    adjlist = Vector{Vector{Tuple{Int,Int}}}()
    rows = rowvals(A)
    for node = 1:size(A,1)
        inds = nzrange(A,node)
        push!(adjlist,collect(zip(rows[inds],A.nzval[inds])))
        # adjlist[node] = collect(zip(rows[inds],A.nzval[inds]))
    end
    return adjlist
end

function _undirected_edges(A::SparseMatrixCSC;use_diag::Bool = false)
    rowval = rowvals(A)
    n = size(A,2)
    ndiag = nnz(diag(A))
    nedges = (nnz(A) - ndiag)/2
    if use_diag
        nedges += ndiag
    end
    edges = Vector{Tuple{Int,Int,Int}}(undef, Int(nedges))
    curedge = 0
    for j=1:n
        for nzi = nzrange(A,j)
            i = rowval[nzi]
            if j > i || (use_diag && i == j)
                curedge += 1
                edges[curedge] = (i,j,A[i,j])
            end
        end
    end
    @assert curedge == nedges
    return edges
end


"""
    traingle_rewire(P::SparseMatrixCSC,offset::Int = 0,ksteps::Int=10000)

TBW
"""
function traingle_rewire(P::SparseMatrixCSC,ksteps::Int=10000;offset::Int = 0)
    #storing the graph in multiple formats for faster access
    adjlist = sparse_to_adjlist(P)
    edges = _undirected_edges(P)

    tri_hist = 0

    for k = 1:ksteps  
        #sample an edge
        edge = rand(1:lastindex(edges))
        edge_i,edge_j,edge_weight = edges[edge]

        #sample a random pair of nodes 
        u,v = extrema(rand(1:lastindex(adjlist),2))
        while u==v || u in first.(adjlist[v]) #nodes equal
            u,v = extrema(rand(1:lastindex(adjlist),2))
        end

        #remove replace old edge with new edge 
        w = commonNeighbors(adjlist,u,v)
        delta_tri = w-edge_weight
        tri_hist += delta_tri
        
        #update edges 
        edges[edge] = (u,v,w)

        #update adjlist
        #remove old edge 
        ind = findfirst(x->x[1]==edge_j, adjlist[edge_i])
        popat!(adjlist[edge_i],ind)

        ind = findfirst(x->x[1]==edge_i, adjlist[edge_j])
        popat!(adjlist[edge_j],ind)

        #add in new edge 
        push!(adjlist[u],(v,w)) 
        push!(adjlist[v],(u,w)) 
    end
    @show tri_hist
    return edges 
end

gname = getgnames("enron","input/graphs/")[1]
A = loadGraph(gname,"input/graphs/")
P = commonNeighbors(A)
adjlist = sparse_to_adjlist(P)

nnodes = size(A,1)
@time new_edges = traingle_rewire(P,nnz(P));

c = 0
for t in triangles(A)
    c += 1
end
c*6

#build new matrix and look at some stats - degrees, triangles, components
ei,ej,ev = Vector{Int}(),Vector{Int}(),Vector{Int}()
for edge in new_edges
    push!(ei,edge[1])
    push!(ej,edge[2])
    push!(ev,edge[3])
end

X = sparse(ei,ej,ones(lastindex(ei)),nnodes,nnodes)
X = max.(X,X')

d = vec(sum(A;dims=1))
d_new = vec(sum(X;dims=1))

histogram(d)
histogram(d_new)

largest_component(X)[1]

c = 0
for t in triangles(X)
    c+=1 
end
@show c 

findall(P.nzval.==0)



#########################################################
############ JUSTIFYING EPIDEMIC PARAMETERS ############# 
#########################################################

hs = ["commutes-all","mexico", "filtered", #row1
    "cn-moduillinois", "cn-Penn", "cn-modWiscon", #row2...
    "dblp","enron","anon", #row 3
    "cit-HepPh", "slashdot", "flickr", 
    "geometric","study-25-150","cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5.smat"]
    #  "rewired-10000.0-modmexico","rewired-10000.0-cn-moduillinois","er-10000.0-dblp-cc"]
gnames = map(h->h=="" ? "" : getgnames(h,"pipeline/graphs/")[1], hs) 
titles = ["US\nCommutes", "Mexico City\nTrace", "Filtered\nUS Flows",
    "Sparsified\nCollege-Illinois", "Sparsified\nCollege-Penn", "Sparsified\nCollege-Wisc.",    
    "Collaboration", "Email", "Facebook\nInteractions",
    "Citation", "Slashdot","Flickr",    
    "Local\nGeometric", "Geometric\nCommunities", "Random Walk\nCommunities"] 

### type 1 plot (full-ish parameters for a single network)

function add_tickmarks(f,xmirror=false,ymirror=false)
    plot!(f,xmirror=xmirror,ymirror=ymirror,
        tick_direction=:out,
        xticks = (range(1.5,24.5,25),["" for i=1:25]),xrotation=90,
        top_margin=2Measures.mm,
        yticks = (range(1.5,15.5,16),["" for i=1:16]),
        colorbar = false,framestyle=true,
        xlims=(1.05,24.95),ylims=(1.05,15.95))
    return f
end

function make_params_plot(gname::String,betas::Vector)
    @assert(lastindex(betas)==15)

    ps = qpercent_contours(gname,betas=betas,colorbar=false,add_labels=false)
    ps = vcat(ps[2]...)   

    for (ind,title) in enumerate(betas)
        f = ps[ind]
    
        ymin,ymax = ylims(f)
        plot!(f,framestyle=:none,
            xlabel="",ylabel="",title="",
            xguidefontsize=12,
            right_margin=0Measures.mm,
            colorbar=false)
        annotate!(f,[(mean(xlims(f)),ymin+0.8*(ymax-ymin),text("β = $(betas[ind])", 16, :center, :center,RGB(0.0,0.0,0.0)))])
    end

    figs = deepcopy(ps)

    for (i,f) in enumerate(figs)
        xmirror,ymirror = false,false
        if i in [1;2;3]
            xmirror = true
        end
        if i in [3;6;9;12;15]
            ymirror = true
        end
        f = add_tickmarks(f,xmirror,ymirror)
        
        if i in 4:12
            plot!(f,xticks=false)
        end
    
        if i in 2:3:15
            plot!(f,yticks=false)
        end
        figs[i] = f
    end
    
    newf = Plots.plot(figs...,
        layout=grid(5,3,heights=[0.2, 0.2, 0.2, 0.2, 0.2],widths=[0.33,0.33,0.33]/(0.33+0.33+0.33)),
        size=(1000,1200),
        # margins=0Measures.mm,
        bottom_margin=-2*Measures.mm,
        right_margin=-1*Measures.mm,
        left_margin=-2*Measures.mm,
        top_margin=-0.5*Measures.mm,
        link=:y)
    

end


gname = getgnames("wisc","input/graphs/")[1]
betas = get_betas(gname)

ps = qpercent_contours(gname,betas=betas,colorbar=false,add_labels=false)
ps = vcat(ps[2]...)

#for tinkering and making adjustments without having to remake figs from scratch
figs = deepcopy(ps[1:15])

#remove text, and relabel
for (ind,title) in enumerate(figs)
    f = figs[ind]

    ymin,ymax = ylims(f)
    plot!(f,framestyle=:none,
        xlabel="",ylabel="",title="",
        xguidefontsize=12,
        right_margin=0Measures.mm,
        colorbar=false)
    annotate!(f,[(mean(xlims(f)),ymin+0.8*(ymax-ymin),text("β = $(betas[ind])", 16, :center, :center,RGB(0.0,0.0,0.0)))])
end

function add_tickmarks(f,xmirror=false,ymirror=false)
    plot!(f,xmirror=xmirror,ymirror=ymirror,
        tick_direction=:out,
        xticks = (range(1.5,24.5,25),["" for i=1:25]),xrotation=90,
        top_margin=2Measures.mm,
        yticks = (range(1.5,15.5,16),["" for i=1:16]),
        colorbar = false,framestyle=true,
        xlims=(1.05,24.95),ylims=(1.05,15.95))
    return f
end


#modifying exterior plots 
for (i,f) in enumerate(figs)
    xmirror,ymirror = false,false
    if i in [1;2;3]
        xmirror = true
    end
    if i in [3;6;9;12;15]
        ymirror = true
    end
    f = add_tickmarks(f,xmirror,ymirror)
    
    if i in 4:12
        plot!(f,xticks=false)
    end

    if i in 2:3:15
        plot!(f,yticks=false)
    end
    figs[i] = f
end
# plot!(figs[2],title=gname[1:end-5],titlefontsize=24)

#plotting
newf = Plots.plot(figs...,
        layout=grid(5,3,heights=[0.2, 0.2, 0.2, 0.2, 0.2],widths=[0.33,0.33,0.33]/(0.33+0.33+0.33)),
        size=(1000,1200),
        # margins=0Measures.mm,
        bottom_margin=-2*Measures.mm,
        right_margin=-1*Measures.mm,
        left_margin=-2*Measures.mm,
        top_margin=-0.5*Measures.mm,
        link=:y)
plot!(newf,dpi=300)
Plots.savefig(newf,"code/paper-figs/heatmaps/uplot-full-params-$(gname[1:end-5]).png")


### doing this more systematically

"""
    graphSize(gname)

returns number of nodes in giant component
"""
function graphSize(gname)
    A = loadGraph(joinpath("pipeline/graphs/",gname))
    return size(A,1)
end

#load in raw data, aggregate, and store for each network then generate plots
"""
    load_raw_data(gname::String)
given gname, load in uniform diffusion data for rewired variants
"""
function load_raw_plotting_data(gname::String)::Dict{Float64,Matrix}
    gnames = getgnames(gname,"pipeline/graphs/")
    filter!(x->!occursin("triangle-rewired-",x),gnames)

    nnodes = graphSize(gname)

    betas = get_betas(gname)
    result = Dict{Float64,Matrix}()
    @show(gname)
    @showprogress for beta in betas
        data,rewiring_fractions1,rewiring_fractions2 = load_double_rewiring_data(gnames[1],
                        beta=beta,gamma=0.05,dtype="tinfs",method="seir")
        data = aggregate_diffusion_data(hcat(map(x->last.(x),data)...))
        result[beta] = deepcopy(data./nnodes)
    end
    return result
end


#make plot from a single graph 
function full_epi_param_plot(gname::String,data::Dict{String,Dict})
    gdata = data[gname]
    betas = keys(gdata)

    #generate contour plots 
    figs = Dict()
    for beta in betas
        beta_data = gdata[beta]
        f = Plots.contourf(1:size(beta_data,2),1:size(beta_data,1),(x,y)->beta_data[y,x],
                c=:heat,clims=(0,1),
                levels=0.0:0.2:1.0)
        figs[beta] = deepcopy(f)
    end
    return figs 
end

#remove labels, colorbar, and annotate plots
function strip_info(gname::String,figs_data::Dict)
    fig_data = figs_data[gname]
    betas = sort(collect(keys(fig_data)))

    ps = Vector()
    for beta in betas
        f = deepcopy(fig_data[beta])
        #remove auxilary info
        ymin,ymax = ylims(f)
        plot!(f,framestyle=:none,
            xlabel="",ylabel="",title="",
            xguidefontsize=12,
            right_margin=0Measures.mm,
            colorbar=false)
        annotate!(f,[(mean(xlims(f)),ymin+0.8*(ymax-ymin),text("β = $beta", 16, :center, :center,RGB(0.0,0.0,0.0)))])
        push!(ps,deepcopy(f))
    end
    return ps
end

#adding tick marks.. need to choose betas before executing these functions 
function add_tickmarks(f,xmirror=false,ymirror=false)
    xbounds,ybounds = xlims(f),ylims(f)
    xmin,xmax = xbounds
    xlen = round(Int,xmax-xmin)+1
    ymin,ymax = ybounds 
    ylen = 16 #16 values of q 

    plot!(f,xmirror=xmirror,ymirror=ymirror,
        tick_direction=:out,
        xticks = (range(xmin+0.5,xmax-0.5,xlen),["" for i=1:xlen]),xrotation=90,
        top_margin=2Measures.mm,
        yticks = (range(ymin+0.5,ymax-0.5,ylen),["" for i=1:ylen]),
        colorbar = false,framestyle=true,
        xlims=(xmin+5e-2,xmax-5e-2),ylims=(1.05,15.95))
    return f
end


function modify_exterior_plots!(ps::Vector,size::Tuple=(7,3))
    #modifying exterior plots 
    for (i,f) in enumerate(ps)
        xmirror,ymirror = false,false
        if i in [1;2;3]
            xmirror = true
        end
        if i in [3;6;9;12;15;18]
            ymirror = true
        end
        f = add_tickmarks(f,xmirror,ymirror)
        
        if i in 4:15
            plot!(f,xticks=false)
        end

        if i in 2:3:18
            plot!(f,yticks=false)
        end
        #fix last plot
        # if i==19
             
        # end

        ps[i] = f
    end
end

function layout_plots(ps::Vector)
    newf = Plots.plot(ps...,
        layout=grid(6,3,heights=[1/6 for i=1:6],widths=[0.33,0.33,0.33]/(0.33+0.33+0.33)),
        size=(1000,1200),
        # margins=0Measures.mm,
        bottom_margin=-2*Measures.mm,
        right_margin=-1*Measures.mm,
        left_margin=-2*Measures.mm,
        top_margin=-0.5*Measures.mm,
        link=:y)
    return newf
end

#putting it all together
hs = ["commutes-all","mexico", "filtered", #row1
    "cn-moduillinois", "cn-Penn", "cn-modWiscon", #row2...
    "dblp","enron","anon", #row 3
    "cit-HepPh", "slashdot", "flickr", 
    "geometric","study-25-150","cl-lfr-100000-3.00-2.00-0.15-1-2000-5-500-17-connect-graph-invdegdepth-8000-0.5-0.0-100-5.smat"]
    #  "rewired-10000.0-modmexico","rewired-10000.0-cn-moduillinois","er-10000.0-dblp-cc"]
gnames = map(h->h=="" ? "" : getgnames(h,"pipeline/graphs/")[1], hs) 

# hs = ["penn","uill","wisc"]
# gnames = map(h->h=="" ? "" : getgnames(h,"pipeline/graphs/")[1], hs) 

# gnames = ["study-11-2023-0-noweak.smat",
#   "study-11-2022-1.smat",
#   "study-11-2022-10.smat",
#   "study-11-2022-20.smat",
#   "study-11-2022-30.smat",
#   "study-11-2022-40.smat",
#   "study-11-2022-50.smat"
#   ]

#load in raw aggreagted data for plotting
data = Dict{String,Dict}()
for gname in gnames
    data[gname] = load_raw_plotting_data(gname)
end

#generate plots 
figs_data = Dict{String,Dict}()
for gname in gnames
    figs_data[gname] = full_epi_param_plot(gname,data)
end

#single graph as an example 
gname = gnames[15]

ps = strip_info(gname,figs_data)
ps = deepcopy(ps)
modify_exterior_plots!(ps[9:end]) #implicitly picking values of beta
newf = layout_plots(ps[9:end]) #actually lay them out 

plot!(newf,dpi=1000)
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/heatmaps/full-params-$(gname[1:end-5]).png")
Plots.savefig(newf,"/p/mnt/scratch/network-epi/code/paper-figs/heatmaps/full-params-$(gname[1:end-5]).pdf")


#just to check every rock..
gname = "study-25-150"
gname = getgnames(gname,"input/graphs/")[1]
ps = qpercent_contours(gname,betas=[2e-2],colorbar=false)
ps = vcat(ps...)
ps[2]

plot!.(ps,framestyle=:none,colorbar=false)
plot(ps[20:29]...,layout=(5,2),
    margins=-1Measures.mm,
    size=(1200,1000))



#########################################################################################
################################# MODIFIED STUDY UPLOTS #################################
#########################################################################################

gs = [ "study-11-2023-0-noweak.smat",
    "study-11-2022-1.smat",
    # "study-11-2022-10.smat",
    # "study-11-2022-20.smat",
    # "study-11-2022-30.smat",
    # "study-11-2022-40.smat",
    # "study-11-2022-45.smat",
    "study-11-2022-50.smat",
]

qcaps = [0,25,50,75,100,125,150,200,250,300,400,500,1000,2000,3000,4000,5000]

#betas for base graph 
bs = unique(vcat(1e-3:1e-3:1e-2,2e-2:1e-2:1e-1))

Abase = loadGraph(gs[1],"input/graphs/")
davg_base = nnz(Abase)/lastindex(Abase,1)

tmp_figs = []
@showprogress for gname in gs 
    #get new beta 
    Atmp = loadGraph(gname,"input/graphs/")
    davg_tmp = nnz(Atmp)/lastindex(Atmp,1)
    

    ps = qpercent_contours(gname,betas=bs.*(davg_base/davg_tmp),colorbar=false,qpercents=vcat(1:lastindex(qcaps)),
            diffusion_type="uniform/scratch")
    ps = vcat(ps[2]...)

    #add verical lines for base network 
    # map(x->plot!(x,[12.5;12.5],[1;17.0],c=:black,linewidth=3),ps)
    # map(x->plot!(x,[13.5;13.5],[1;17.0],c=:black,linewidth=3,leg=false),ps)

    # map(x->plot!(x,colorbar=false,ylims=(1,17),xlims=(0.95,25.05)),ps)
    
    push!(tmp_figs,deepcopy.(ps)...)
end


ind = 10
tt = vcat(tmp_figs[ind:ind+5],tmp_figs[ind+19:ind+5+19],tmp_figs[ind+38:ind+5+38])
f = tt[1]
plot!(f,yicks=(1:17,string.(qcaps)),xticks=false)

map(x->plot!(x,yaxis=false),tt[2:end])

tt[2]

newf = plot(deepcopy(tt)...,layout=(3,6),colorbar=false,
    size=(1800,1000),
    bottom_margin=-5Measures.mm,
    top_margin=-2Measures.mm,
    right_margin=-1.5Measures.mm,
    left_margin=-10Measures.mm)





plot(tmp_figs[ind],tmp_figs[ind+19],tmp_figs[ind+38],layout=(1,3),
    xticks=false,yticks=false,size=(1000,300),
    right_margin=-1Measures.mm,
    left_margin=-1Measures.mm,
    bottom_margin=-1Measures.mm,
    top_margin=-1Measures.mm)

plot(tmp_figs...,layout=(2,4),xaxis=false,yaxis=false,
        left_margin=-8Measures.mm,
        right_margin=-2Measures.mm,
        bottom_margin=-4Measures.mm,
        top_margin=-2Measures.mm,
        size=(1600,800))


####
gnames =[
    "internal-shuffled-cl-louvain-1-1-study-25-1.smat",
    "internal-shuffled-cl-louvain-2-1-study-25-1.smat",
    "internal-shuffled-cl-louvain-3-1-study-25-1.smat"
]


g = gnames[1]
rs = qpercent_contours(g,betas=get_betas(g),colorbar=false,add_labels=false)[2]

g = gnames[2]
ts = qpercent_contours(g,betas=get_betas(g),colorbar=false,add_labels=false)[2]

g = gnames[3]
us = qpercent_contours(g,betas=get_betas(g),colorbar=false,add_labels=false)[2]

fs = []
push!(fs,deepcopy(rs[15]),deepcopy(ts[15]),deepcopy(us[15]))

plot(fs...,colorbar=false,xticks=false,yticks=false,layout=(1,3))
b = get_betas(gname)[15]

function add_tickmarks(f,xmirror=false,ymirror=false)
    plot!(f,xmirror=xmirror,ymirror=ymirror,
        tick_direction=:out,
        xticks = (range(1.5,24.5,25),["" for i=1:25]),xrotation=90,
        top_margin=2Measures.mm,
        yticks = (range(1.5,15.5,16),["" for i=1:16]),
        colorbar = false,framestyle=true,
        xlims=(1.05,24.95),ylims=(1.05,15.95))
    return f
end

fs[1] = add_tickmarks(fs[1])
fs[2] = add_tickmarks(fs[2])
fs[3] = add_tickmarks(fs[3])

plot!(fs[1],colorbar=false,
    # title="internal-shuffled-cl-louvain-1-1-$(gname[1:end-5])\nβ=$b",
    topmargin=3Measures.mm,
    dpi=300)
plot!(fs[2],colorbar=false,yticks=false,
    # title="internal-shuffled-cl-louvain-2-1-$(gname[1:end-5])\nβ=$b",
    topmargin=3Measures.mm,
    dpi=300)
plot!(fs[3],colorbar=false,ymirror=true,
    # title="internal-shuffled-cl-louvain-3-1-$(gname[1:end-5])\nβ=$b",
    topmargin=3Measures.mm,
    dpi=300)

Plots.savefig(fs[1],"/p/mnt/scratch/network-epi/code/paper-figs/example-figs/meso-$(gnames[1][1:end-5]).png")
Plots.savefig(fs[2],"/p/mnt/scratch/network-epi/code/paper-figs/example-figs/meso-$(gnames[2][1:end-5]).png")
Plots.savefig(fs[3],"/p/mnt/scratch/network-epi/code/paper-figs/example-figs/meso-$(gnames[3][1:end-5]).png")