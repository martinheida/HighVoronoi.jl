
struct ThreadsafeProgressMeter{RWL<:Union{BusyFIFOLock,Nothing}}
    p::Progress
    silence::Bool
    lock::RWL
end
function ThreadsafeProgressMeter(tpm::TPM) where {TPM<:ThreadsafeProgressMeter}
    return tpm
end

function ThreadsafeProgressMeter(n::Int,silence::Bool,intro) 
    return ThreadsafeProgressMeter(Progress(n,intro),silence,nothing)
end

function ThreadsafeProgressMeter(n::Int,silence::Bool,intro,::MultiThread) 
    return ThreadsafeProgressMeter(Progress(n,intro),silence,BusyFIFOLock())
end
 
#import Progress
@inline function next!(tpm::TPM) where {TPM<:ThreadsafeProgressMeter}
    lock(tpm.lock)
    (!tpm.silence) && (ProgressMeter.next!(tpm.p))
    unlock(tpm.lock)
end
#ProgressMeter.lock_if_threading(f::Function, p::ProgressMeter.Progress) = f()


###################################################################################################

## provides a bunch of tools to display progress of voronoi algorithms properly

###################################################################################################

### come offset functions

const voronoi_offset    = 4
const integral_offset   = 4
const sys_refine_offset = 4
const BC_offset         = 4

#### Main functions

function vp_line()
    print("\n")
end

#=function vp_column(i)
    print("\u1b[0E\u1b[$(i)C")    
end

function vp_delete_from_here()
    
end

function vp_delete_line_content()
    print("\u1b[2K") # delete entire line
end=#

function vp_line_up()
    print("\u1b[2K\u1b[1A\u1b[200D")
end

#=
function vp_line_up(K)
    for i in 1:K
        print("\u1b[2K\u1b[1A\u1b[200D")
    end
end

function vp_blocks(content,offsets)
    for i in 1:(length(content))
        print("\u1b[0E\u1b[$(offset[i])C")    
        print(content[i])
    end
end
=#

function vp_print(o1::Int,c;crayon=nothing)
#    if typeof(crayon)==Nothing
#        print("\u1b[0E\u1b[$(o1)C")    
        print("\u1b[200D\u1b[$(o1)C")    
        print(c)
#    else
#=        print(crayon)
        print("\u1b[0E\u1b[$(o1)C")    
        print(c)
        print(Crayon(reset=true))=#
#    end
end

function vp_print(o1::Int,c1,o2::Int,c2;crayon=nothing)
#    if typeof(crayon)==Nothing
        print("\u1b[0E\u1b[$(o1)C")    
        print(c1)
        print("\u1b[0E\u1b[$(o2)C")    
        print(c2)
#    else
#=        print(crayon)
        print("\u1b[0E\u1b[$(o1)C")    
        print(c1)
        print("\u1b[0E\u1b[$(o2)C")    
        print(c2)
        print(Crayon(reset=true))=#
#    end
end
