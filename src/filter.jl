struct Filter{M <: AbstractMesh, F<:Function, II<:AbstractVector{Int64}}
    mesh::M
    iterate::II
    affected::BitVector
    condition::F
    lm::Int64
    lm_internal::Int64
    function Filter(m::M, i::II,f::F) where {M <: AbstractMesh, F<:Function, II<:AbstractVector{Int64}}
        _affected = falses(length(m))
        _affected[i] .= true
        new{M,F,II}(m, i, _affected, f, length(m),0)#internal_length(m))
    end
end

@inline function mark_sig(_filter,sig) 
    for s in sig
        (s>_filter.lm || s==0) && continue
        _filter.affected[s] = true
    end
end


function filter_data( condition, mesh,_depsig,_depr, affected) 
    depsig = StaticBool(_depsig)
    depr = StaticBool(_depr)
    parameter(sig,r,dep_sig::StaticTrue,dep_r::StaticTrue) = (sig,r)
    parameter(sig,r,dep_sig::StaticTrue,dep_r::StaticFalse) = (sig,)
    parameter(sig,r,dep_sig::StaticFalse,dep_r::StaticTrue) = (r,)
    f = (sig,r)->condition(parameter(sig,r,depsig,depr)...)
    getaffected(aff::AbstractVector{Bool},m) = view(1:length(m),aff)
    getaffected(aff,m) = aff
    getaffected(aff::UnitRange{Int},m) = collect(aff)
    myaffected = getaffected(affected,mesh)
    _filter = Filter(mesh,myaffected,condition)
    return myaffected,_filter
end

#=
filtermesh!(mesh,my_affected,_filter)
cleanupfilter!(mesh,internal_index(mesh,i))
mark_sig(_filter,sig)
mark_delete_vertex!(_filter,vdb,sig)

mark_delete_bv!(mesh,edge,bv)
cleanupfilter_bv!(mesh)

boundary_vertices_iterator(mesh)
=#
@inline mark_delete_ref(_filter,sig,ref::Nothing) = nothing
#@inline mark_delete_ref(_filter,sig,ref::Int) = nothing

function mark_delete_ref(_filter,sig,ref::U) where {U<:Union{VertexRef,Int}}
    il = internal_length(_filter.mesh)
    for i in 2:length(sig)
        s = sig[i]
        s>il && continue
        delete_reference(_filter.mesh,s,ref)
    end
end

function filter!( condition, mesh::M,_depsig=StaticBool{true}(),_depr=StaticBool{true}(), _sig_internal=!_depsig; affected = 1:length(mesh),searcher=nothing) where {M<:AbstractMesh} 
    my_affected, _filter = filter_data( condition, mesh,_depsig,_depr, affected)
    range = typeof(searcher)==Nothing ? (1:2) : ((length(mesh)+1):length(mesh)+length(searcher.domain))
    #println("filter start")
    for i in my_affected
        #print("1") 
        activate_cell(searcher,i,range)
        avi = all_vertices_iterator(_filter.mesh,internal_index(_filter.mesh,i),statictrue)
        #print("2")
        for (sig,r) in avi
            #print("a")
            ext_sig = external_sig(_filter.mesh,sig,statictrue)
            #print("b")
            if _sig_internal==statictrue ? !_filter.condition(sig,r) : !_filter.condition(ext_sig,r)
                #print("A")
                mark_sig(_filter,ext_sig)
                #print("B")
                ref = mark_delete_vertex!(_filter.mesh,sig,IterationIndex(avi))
                #print("C")
                mark_delete_ref(_filter,sig,ref)
                #print("D")
            end
            #print("c")
        end
        #print("3")
    end
    #println("filter mitte")
    for i in 1:_filter.lm
        _filter.affected[i] && cleanupfilter!(mesh,internal_index(mesh,i)) # provide second arg in internal representation
    end
    #println("filter ende")
    #return my_affected
end

""" filters verteces of `affected` according to `condition`. Does NOT reduce number of points. This must be done manually"""
function filter_bv!( condition, mesh::M,_depsig=StaticBool{true}(),_depr=StaticBool{true}(); affected = nodes_iterator(mesh) ) where {M<:AbstractMesh}
    myaffected, _filter = filter_data( condition, mesh,_depsig,_depr, affected)
    for (edge,bv) in boundary_vertices_iterator(mesh)
        if !_filter.condition(edge,bv) # if sig consists only of affected cells
            mark_delete_bv!(mesh,edge,bv)
        end
    end
    cleanupfilter_bv!(mesh)
end


