
# Definiere einen abstrakten Datentyp HVDataBase
abstract type HVDataBase{P<:Point, IV<:AbstractVector{Int64}} end

PointType(::Type{HVD}) where {P, IV, HVD<:HVDataBase{P, IV}} = P
VectorType(::Type{HVD}) where {P, IV, HVD<:HVDataBase{P, IV}} = IV
PointType(::HVD) where {P, IV, HVD<:HVDataBase{P, IV}} = P
VectorType(::HVD) where {P, IV, HVD<:HVDataBase{P, IV}} = IV
blocksize(::Type{HVD}) where {P, IV, HVD<:HVDataBase{P, IV}} = begin
    dimension = size(P)[1]
    return (1+2*dimension) * round(Int,lowerbound(dimension,dimension))
end
blocksize(::HVD) where {P, IV, HVD<:HVDataBase{P, IV}} = begin
    dimension = size(P)[1]
    return (1+2*dimension) * round(Int,lowerbound(dimension,dimension))
end

################ Expected Functionality:

##  push!(::HKD, sig=>r) : returns Int64 - Index
##  get(::HKD, index) : returns sig=>r for index
##  haskey(::HKD, sig) : returns whether has key or not.


##############################################################################################################################################################
##############################################################################################################################################################

## HeapDataBase

##############################################################################################################################################################
##############################################################################################################################################################


mutable struct HeapDataBase{P, RI, Q<:QueueHashTable} <: HVDataBase{P, Vector{Int64}}
    data::Vector{Vector{Int64}}
    floatview::Vector{RI}
    keys::Q
    blocksize::Int64
    position::Int64
    buffer_sig::Vector{Vector{Int64}}
    buffer_r::Vector{Vector{Float64}}
    lock::ReadWriteLock
    nthreads::Int64
    set::Dict{Vector{Int64},Int64}

    # Konstruktor für HeapDataBase
    function HeapDataBase(::Type{P},_blocksize::Int64) where {P}
        data = Vector{Vector{Int64}}() # Initialisiere leeres data
        RI = typeof(Base.reinterpret(Float64,Int64[]))
        floatview = Vector{RI}()       # Initialisiere leeres floatview
        queue = QueueHashTable(_blocksize,SingleThread())
        nth = Threads.nthreads()
        #new{P, RI, typeof(queue)}(data, floatview, queue, _blocksize,_blocksize #=position=blocksize??=#,[Int64[] for i in 1:nth],[Float64[] for i in 1:nth],ReadWriteLock(),nth)
        new{P, RI, typeof(queue)}(data, floatview, queue, _blocksize,_blocksize,[Int64[] for i in 1:nth],[Float64[] for i in 1:nth],ReadWriteLock(),nth,Dict{Vector{Int64},Int64}())
    end
end
indexvector(::HDB) where {HDB<:HeapDataBase} = Vector{Int64}()

@inline function push_entry!(hdb::HDB,d::Float64) where {HDB<:HeapDataBase}
    hdb.position +=1 
    if hdb.position>hdb.blocksize
        push!(hdb.data,Vector{Int64}(undef,hdb.blocksize))
        push!(hdb.floatview,Base.reinterpret(Float64,hdb.data[end]))
        hdb.position = 1
    end
    hdb.floatview[end][hdb.position] = d
end

@inline function push_entry!(hdb::HDB,d::Int64) where {HDB<:HeapDataBase}
    hdb.position +=1 
    if hdb.position>hdb.blocksize
        push!(hdb.data,Vector{Int64}(undef,hdb.blocksize))
        push!(hdb.floatview,Base.reinterpret(Float64,hdb.data[end]))
        hdb.position = 1
    end
    hdb.data[end][hdb.position] = d
end

function push!(hdb::HDB,d::AV) where {HDB<:HeapDataBase,R<:Real,AV<:AbstractVector{R}}
    @inline db(::Type{Float64},hd) = hd.floatview[end]
    @inline db(::Type{Int64},hd) = hd.data[end]
    l = length(d)
    if (l>hdb.blocksize-hdb.position)
        for dd in d
            push_entry!(hdb,dd)
        end
    else
        mydb = db(R,hdb)
        view(mydb,(hdb.position+1):(hdb.position+l)) .= d
        hdb.position += l
    end
end

function push!(hdb::HDB,pair::P) where {HDB<:HeapDataBase,P<:Pair}
    sig = pair[1]
    r = pair[2]
    writelock(hdb.lock)
    pushqueue!(hdb.keys,sig)
    push_entry!(hdb,length(sig))
    pos = hdb.position + hdb.blocksize*(length(hdb.data)-1)
    push!(hdb,sig)
    push!(hdb,r)
    writeunlock(hdb.lock)
    if pos==0
        open("protokol.txt","a") do f
            write(f,"hier geht's shon los in push!(HDB) \n")
        end
    end
    return pos
end
@inline haskey(hdb::HDB,key::K) where {HDB<:HeapDataBase,K} = begin
    readlock(hdb.lock)
    r = haskey(hdb.keys,key)
    readunlock(hdb.lock)
    return r
end

function read_from_db(data,buffer_,pos,block,l,blocksize)
    if length(buffer_)<l
        resize!(buffer_,l)
    end
    sig = view(buffer_,1:l)
    l1 = blocksize-pos
    if (l<=l1)
        sig .= view(data[block],(pos+1):(pos+l))
    else
        view(buffer_,1:l1) .= view(data[block],(pos+1):(pos+l1))
        #try
            view(buffer_,(l1+1):l) .= view(data[block+1],1:(l-l1))
        #catch
        #    error("$(length(data[block+1])), $blocksize")
        #end
    end
    return sig
end
function rescale_db(hdb::HDB) where {P,HDB<:HeapDataBase{P}}
    writelock(hdb.lock)
    nth = Threads.nthreads()
    resize(hdb.buffer_sig,nth)
    resize(hdb.buffer_r,nth)
    for i in (hdb.nthreads+1):nth
        hdb.buffer_sig = Int64[]
        hdb.buffer_r = Float64[]
    end
    hdb.nthreads = nth
    writeunlock(hdb.lock)
end
function get_entry(hdb::HDB,index) where {P,HDB<:HeapDataBase{P}}
    id = Threads.threadid()
    id>hdb.nthreads && rescale_db(hdb)
    pos = mod(index,hdb.blocksize)
    block = div(index,hdb.blocksize)+1
    if pos==0
        pos = hdb.blocksize
        block -= 1
    end
    readlock(hdb.lock)
    #if pos==0 || block==0 
    #    open("protokol.txt","a") do f
    #    write(f,"Pos: $pos, Block: $block, Index: $index \n")
    #    end
    #end
    @inbounds l = hdb.data[block][pos]
    if l<0 
        readunlock(hdb.lock)
        return (view(hdb.buffer_sig[id],1:0),zeros(P)) 
    end
    sig = read_from_db(hdb.data,hdb.buffer_sig[id],pos,block,l,hdb.blocksize)
    pos2 = pos + l
    if (pos2>=hdb.blocksize)
        pos2 -= hdb.blocksize
        block += 1
    end
    r = P(read_from_db(hdb.floatview,hdb.buffer_r[id],pos2,block,length(P),hdb.blocksize))
    readunlock(hdb.lock)
    return (sig,r)
end

function Base.delete!(hdb::HDB,index,key) where {P,HDB<:HeapDataBase{P}}
    pos = mod(index,hdb.blocksize)
    block = div(index,hdb.blocksize)+1
    if pos==0
        pos = hdb.blocksize
        block -= 1
    end
    #print("*")
    writelock(hdb.lock)
    hdb.data[block][pos]*=-1
    delete!(hdb.keys,key)
    writeunlock(hdb.lock)
    #print("/")
end 


function push_facet!(hdb::HDB,data::P) where {HDB<:HeapDataBase,P<:Tuple}
    sig = data[1]
    r = data[2]
    u = data[3]
    sort!(sig)
    pos = 0
    scs = copy(sig)
    #if haskey(hdb.set,scs)
    #    println()
    #    println("----- KENN ICH SCHON !!! ")
    #end
    #push!(hdb.set,scs=>1)
    status(41)
    writelock(hdb.lock)
    if !(pushqueue!(hdb.keys,sig)) # works like a haskey but also pushes entry if not known
        status(42)
        status(43)
        push_entry!(hdb,length(sig))
            pos = hdb.position + hdb.blocksize*(length(hdb.data)-1)
            push!(hdb,sig)
            push!(hdb,r)
            push!(hdb,u)
            status(44)
    else
        status(45)
        #println()
        #println("KENN ICH SCHON !!! ")
    end
    writeunlock(hdb.lock)
    return pos
end


function get_facet(hdb::HDB,index) where {P,HDB<:HeapDataBase{P}}
    id = Threads.threadid()
    index==0 && error("severe error: index shall never be '0' !")
    id>hdb.nthreads && rescale_db(hdb)
    pos = mod(index,hdb.blocksize)
    block = div(index,hdb.blocksize)+1
    if pos==0
        pos = hdb.blocksize
        block -= 1
    end
    readlock(hdb.lock)
    #if pos==0 || block==0 
    #    open("protokol.txt","a") do f
    #    write(f,"Pos: $pos, Block: $block, Index: $index \n")
    #    end
    #end
    #@inbounds 
    l = hdb.data[block][pos]
    if l<0 
        readunlock(hdb.lock)
        return (view(hdb.buffer_sig[id],1:0),zeros(P)) 
    end
    sig = read_from_db(hdb.data,hdb.buffer_sig[id],pos,block,l,hdb.blocksize)
    pos2 = pos + l
    if (pos2>=hdb.blocksize)
        pos2 -= hdb.blocksize
        block += 1
    end
    data = read_from_db(hdb.floatview,hdb.buffer_r[id],pos2,block,2*length(P),hdb.blocksize)
    r = P(view(data,1:length(P)))
    u = P(view(data,(length(P)+1):(2*length(P))))
    readunlock(hdb.lock)
    return (sig,r,u)
end

