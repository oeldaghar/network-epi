# fast chunglu see https://scholar.harvard.edu/files/joelmiller/files/fast_chung_lu_generation.pdf

using SparseArrays
using MatrixNetworks

function fast_chung_lu_generation(wvec)
    @assert(all(wvec.>0),"node weights must be positive")
    @assert(wvec[i]>=wvec[i+1] for i=1:length(wvec)-1)
    ei,ej = zeros(Int),zeros(Int)
    s = sum(wvec)
    nnodes = length(wvec)
    for u = 1:nnodes-1
        v = u+1
        wu=wvec[u]
        p = min(wu*wvec[v]/s,1)
        while v<=nnodes && p>0
            if p!=1
                r=rand()
                v=v+floor(Int,log(r)/log(1-p))
            end
            if v<=nnodes
                q = min(wu*wvec[v]/s,1)
                r = rand()
                if r < q/p
                    push!(ei,u)
                    push!(ej,v)
                end
                p = q
                v+=1
            end
        end
    end
    B = sparse(ei,ej,ones(length(ei)),nnodes,nnodes)
    return max.(B,B')
end

#need to adjust this since our weights are not sorted in decending order but need to be for
# this algorithm to work.
function fast_chung_lu_connect_components(A::SparseMatrixCSC)
    comps = scomponents(A)
    if comps.number != 1 #more than one component so do chung lu Model
        ei,ej,ev = findnz(A)
        wvec = vec(sum(A;dims=2))
        nodeperm = sortperm(wvec,rev=true)
        invnodeperm = sortperm(nodeperm)

        wvec = wvec[nodeperm] #sorted in descending order. necessary for algorithm

        s = sum(wvec)
        nnodes = length(wvec)
        cmap = collect(1:comps.number)
        for i = 1:nnodes-1
            u = nodeperm[i]
            cu = comps.map[u] #original component of node u
            v = u+1
            wu=wvec[u]
            p = min(wu*wvec[v]/s,1)
            while v<=nnodes && p>0
                if p!=1
                    r=rand()
                    v=v+floor(Int,log(r)/log(1-p))
                end
                if v<=nnodes
                    cv = comps.map[v]
                    q = min(wu*wvec[v]/s,1)
                    r = rand()
                    if r < q/p && cmap[cu]!=cmap[cv] #u and v in diff components
                        push!(ei,u)
                        push!(ej,v)
                        push!(ei,v)
                        push!(ej,u)
                        #do component update
                        newcv = cmap[cv]
                        newcu = cmap[cu]

                        c = min(newcu,newcv)
                        @inbounds for i=1:length(cmap)
                            if cmap[i]==newcu || cmap[i]==newcv
                                cmap[i]=c
                            end
                        end

                        # cmap[cmap.==cv].=cmap[cu]
                        # println("$cu   $(cmap[comps.map[v]])")
                        # comps.map[(comps.map).==comps.map[u]].=comps.map[v]
                        comps.number-=1
                    end
                    p = q
                    v+=1
                end
            end
        end #end edge insertions
        A = sparse(ei,ej,ones(length(ei)),size(A)...)
    end
    return A#,cmap
end

#TODO write the above fcn to keep node labels intact
# function fast_chung_lu_connect_components1(A::SparseMatrixCSC)
#     comps = scomponents(A)
#     cmap = deepcopy(comps.map)
#     if comps.number != 1 #more than one component so do chung lu Model
#         ei,ej,ev = findnz(A)
#         wvec = vec(sum(A;dims=2))
#
#         s = sum(wvec)
#         nnodes = length(wvec)
#
#         for u = 1:nnodes-1
#             cu = comps.map[u]
#             v = u+1
#             wu=wvec[u]
#             p = min(wu*wvec[v]/s,1)
#             while v<=nnodes && p>0
#                 if p!=1
#                     r=rand()
#                     v=v+floor(Int,log(r)/log(1-p))
#                 end
#                 if v<=nnodes
#                     cv = comps.map[v]
#                     q = min(wu*wvec[v]/s,1)
#                     r = rand()
#                     if r < q/p && cmap[cu]!=cmap[cv] #u and v in diff components
#                         push!(ei,u)
#                         push!(ej,v)
#                         push!(ei,v)
#                         push!(ej,u)
#                         #do component update
#                         c = min(cmap[cv],cmap[cu])
#                         for i=1:length(cmap)
#                             if cmap[i]==cv || cmap[i]==cu
#                                 cmap[i]=c
#                             end
#                         end
#
#                         # cmap[cmap.==cv].=cmap[cu]
#                         # println("$cu   $(cmap[comps.map[v]])")
#                         # comps.map[(comps.map).==comps.map[u]].=comps.map[v]
#                         comps.number-=1
#                     end
#                     p = q
#                     v+=1
#                 end
#             end
#         end #end edge insertions
#         A = sparse(ei,ej,ones(length(ei)),size(A)...)
#     end
#     return A,cmap
# end

function fast_chung_lu_connect_components1(A::SparseMatrixCSC)
    comps = scomponents(A)
    if comps.number != 1 #more than one component so do chung lu Model
        ei,ej,ev = findnz(A)
        wvec = vec(sum(A;dims=2))
        nodeperm = sortperm(wvec,rev=true)
        invnodeperm = sortperm(nodeperm)

        wvec = wvec[nodeperm] #sorted in descending order. necessary for algorithm

        s = sum(wvec)
        nnodes = length(wvec)
        cmap = collect(1:comps.number)
        for i = 1:nnodes-1
            u = nodeperm[i]
            cu = comps.map[u] #original component of node u
            v = u+1
            wu=wvec[u]
            p = min(wu*wvec[v]/s,1)
            while v<=nnodes && p>0
                if p!=1
                    r=rand()
                    v=v+floor(Int,log(r)/log(1-p))
                end
                if v<=nnodes
                    cv = comps.map[v]
                    q = min(wu*wvec[v]/s,1)
                    r = rand()
                    if r < q/p && cu!=cv #u and v in diff components
                        push!(ei,u)
                        push!(ej,v)
                        push!(ei,v)
                        push!(ej,u)
                        #do component update
                        # newcv = cmap[cv]
                        # newcu = cmap[cu]

                        # c = min(newcu,newcv)
                        # @inbounds for i=1:length(cmap)
                        #     if cmap[i]==newcu || cmap[i]==newcv
                        #         cmap[i]=c
                        #     end
                        # end

                        # cmap[cmap.==cv].=cmap[cu]
                        # println("$cu   $(cmap[comps.map[v]])")
                        # comps.map[(comps.map).==comps.map[u]].=comps.map[v]
                        # comps.number-=1
                    end
                    p = q
                    v+=1
                end
            end
        end #end edge insertions
        A = sparse(ei,ej,ones(length(ei)),size(A)...)
    end
    return A#,cmap
end
