

###############################################################################################################################

## AbstractFile

###############################################################################################################################
#=
using StaticArrays
struct StaticBool{S}
    function StaticBool{S}() where {S}
        new{S::Bool}()
    end
end
Base.@pure StaticBool(S::Bool) = StaticBool{S}()
Base.@pure StaticBool(S::StaticBool{true}) = StaticBool{true}()
Base.@pure StaticBool(S::StaticBool{false}) = StaticBool{false}()
const StaticTrue = StaticBool{true}
const StaticFalse = StaticBool{false}
const statictrue = StaticTrue()
const staticfalse = StaticFalse()
#Base.@pure StaticBool(S) = StaticBool{false}()
@inline Base.:(==)(x::StaticBool{true}, y::Bool) = y == true
@inline Base.:(==)(x::StaticBool{false}, y::Bool) = y == false
@inline @generated Base.:(==)(y::Bool, x::SB ) where SB<:StaticBool = :(x==y)

@inline  Base.:(==)(x::StaticBool{false}, y::StaticBool{false}) = true
@inline  Base.:(==)(x::StaticBool{false}, y::StaticBool{true}) = false
@inline  Base.:(==)(x::StaticBool{true}, y::StaticBool{false}) = false
@inline  Base.:(==)(x::StaticBool{true}, y::StaticBool{true}) = true
@inline  Base.:(!)(x::StaticBool{true}) = staticfalse
@inline  Base.:(!)(x::StaticBool{false}) = statictrue
@inline Base.Bool(x::StaticBool{A}) where A = A



mutable struct ShortVector{T} <: AbstractVector{T}
    data::T
end
ShortVector{T}() where {T<:Real} = ShortVector{T}(0)
ShortVector{T}() where {T} = ShortVector{T}(T())

Base.length(::ShortVector{T}) where {T} = 1
Base.size(::ShortVector{T}) where {T} = (1,)
Base.getindex(v::ShortVector{T}, i::Int) where {T} = (i == 1) ? v.data : throw(BoundsError(v, i))
Base.setindex!(v::ShortVector{T}, value::T, i::Int) where {T} = (i == 1) ? (v.data = value) : throw(BoundsError(v, i))
Base.iterate(v::ShortVector{T}, state=1) where {T} = state == 1 ? (v.data, 2) : nothing
Base.show(io::IO, v::ShortVector{T}) where {T} = print(io, "ShortVector(", v.data, ")")
=#
abstract type AbstractFile end
# File types built upon HVFile
mutable struct MutableFile{HVF<:Union{AbstractFile,IO}} <: AbstractFile
    file::HVF
end
struct PermanentFile{HVF<:AbstractFile} <: AbstractFile
    file::HVF
end

@inline Base.seek(file::AF,position) where {AF<:AbstractFile} = seek(file.file,position) 
@inline Base.seekend(file::AF) where {AF<:AbstractFile} = seekend(file.file) 
@inline Base.position(file::AF) where {AF<:AbstractFile} = position(file.file) 
@inline Base.read(file::AF,type) where {AF<:AbstractFile} = read(file.file,type) 
@inline Base.write(file::AF,data) where {AF<:AbstractFile} = write(file.file,data) 
Base.close(file::AF) where {AF<:AbstractFile} = close(file.file)
function allocate(file::HVF,size::Int64,free::Bool;data = nothing) where {HVF<:Union{AbstractFile,IO}}
    wfirst(::Nothing) = nothing
    wfirst(a) = write(file,a)
    if free
        seekend(file)
        pos = position(file)
        wfirst(data)
        write(file,zeros(Int8,size))
        return pos
    else
        allocate(file,size,data=data)
    end    
end

@inline function write_Int64_position(file::AF,position,i,val::Int64) where {AF<:Union{AbstractFile,IO}}
    seek(file,position + (i-1)*sizeof(Int64))
    write(file,val)
end
@inline function write_Int64_position(file::AF,position,i,vals::AVI) where {AF<:Union{AbstractFile,IO},AVI<:AbstractVector{Int64}}
    seek(file,position + (i-1)*sizeof(Int64))
    for val in vals
        write(file,val)
    end
end

fdfa_size(::Type{T}) where {T<:Real} = sizeof(T)



################################################################################################################

## Blocks

################################################################################################################

struct Block{T}
    data::Vector{T}
    parameters::MVector{2,Int64}
    Block{T}() where {T} = new{T}(Vector{T}(undef,0),MVector{2,Int64}([0,0]))
    function Block{T}(array) where {T}
        b = Block{T}()
        readblock(array.file,array.first,b)
    end
end
@inline function Base.empty!(b::Block)
    resize!(b.data,0)
    getfield(b,:parameters).=0
end
@inline Base.getproperty(cd::Block, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::Block, ::Val{:length}) =  :(getfield(cd,:parameters)[1])
@inline @generated dyncast_get(cd::Block, ::Val{:next}) =  :(getfield(cd,:parameters)[2])
@inline @generated dyncast_get(cd::Block, ::Val{:data}) =  :(getfield(cd,:data))
@inline Base.setproperty!(cd::Block, prop::Symbol, val) = dyncast_set(cd,Val(prop),val)
@inline @generated dyncast_set(cd::Block, ::Val{:length}, val) = :(getfield(cd,:parameters)[1]=val)
@inline @generated dyncast_set(cd::Block, ::Val{:next}, val) = :(getfield(cd,:parameters)[2]=val)
@inline Base.getindex(b::Block{T},index::Int64) where T = b.data[index]

function readblock(file,position,data::Block{T},traceless=false) where T
    pos = traceless ? 0 : Base.position(file)
    if position>0 
        seek(file,position)
        data.length=read(file,Int64)
        data.next = read(file,Int64)
    else 
        data.length = 0
        data.next = 0
    end
    resize!(data.data,data.length)
    for k in 1:data.length
        data.data[k] = read(file,T)
    end
    traceless && seek(file,pos)
    return data
end
@inline readblock(file,position,::Type{T}) where T = readblock(file,position,Block{T}())

function readblock(file,data::Block{T},length::Int64) where T
    data.length=length
    data.next = read(file,Int64)
    resize!(data.data,data.length)
    for k in 1:data.length
        data.data[k] = read(file,T)
    end
    return data
end


################################################################################################################

## FileArrayIterator

################################################################################################################


struct FileArrayIterator{T,A}
    block::Block{T}
    data::MVector{2,Int64}
    array::A
    function FileArrayIterator{T}(array::A) where {T,A}
        fai = new{T,A}(Block{T}(),MVector{2,Int64}(zeros(Int64,2)),array)
        readblock(array.file,array.first,fai.block)
        return fai
    end
end

function increase(fai::FileArrayIterator{T}) where T
    fai.inblock = fai.inblock + 1
    fai.index = fai.index + 1
    if fai.index>fai.array.offset && fai.block.length>0
        empty!(fai.block)
        fai.inblock = 1
    end
    if fai.inblock>fai.block.length && fai.block.length>0
            readblock(fai.array.file,fai.block.next,fai.block)
            fai.inblock = 1
    end
end

# may only be called if fai.valid==true
function get_current(fai::FileArrayIterator{T}) where T
    if fai.inblock>fai.block.length
        return fai.array.buffer[fai.inblock]
    else 
        return getfield(fai,:block).data[getfield(fai,:data)[1]]
    end
end

@inline Base.getproperty(cd::FileArrayIterator, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::FileArrayIterator, ::Val{:inblock}) =  :(getfield(cd,:data)[1])
@inline @generated dyncast_get(cd::FileArrayIterator, ::Val{:index}) =  :(getfield(cd,:data)[2])
@inline @generated dyncast_get(cd::FileArrayIterator, ::Val{:current}) =  :(get_current(cd))
@inline @generated dyncast_get(fai::FileArrayIterator, ::Val{:valid}) =  :(fai.inblock<=fai.block.length || fai.inblock<=fai.array.bufferposition)
@inline @generated dyncast_get(cd::FileArrayIterator, ::Val{S}) where S = :( getfield(cd, S))
@inline Base.setproperty!(cd::FileArrayIterator, prop::Symbol, val) = dyncast_set(cd,Val(prop),val)
@inline @generated dyncast_set(cd::FileArrayIterator, ::Val{:inblock}, val) = :(getfield(cd,:data)[1]=val)
@inline @generated dyncast_set(cd::FileArrayIterator, ::Val{:index}, val) = :(getfield(cd,:data)[2]=val)

################################################################################################################

## Block-Vectors

################################################################################################################

abstract type AbstractBlockVector{T,F} end

function provide_block(lb::ABV,index::Int,block = Block{T}()) where {T,F,ABV<:AbstractBlockVector{T,F}}
    len,current = find_block(lb,index)
    return readblock(lb.file,block,len), index-current
end

function moveto(lb::ABV,index::Int) where {T,F,ABV<:AbstractBlockVector{T,F}}
    len,current = find_block(lb,index)
    read(lb.file,Int64)
    for i in 1:(index-current-1)
        read(lb.file,T)
    end
end

@inline function Base.getindex(lb::ABV,index::Int64) where {T,F,ABV<:AbstractBlockVector{T,F}}
    moveto(lb,index)
    return read(lb.file,T)
end

@inline function Base.setindex!(lb::ABV,val::T,index::Int64) where {T,F,ABV<:AbstractBlockVector{T,F}}
    moveto(lb,index)
    return write(lb.file,val)
end

################################################################################################################

## Linear Block-Vectors

################################################################################################################

struct LinearBlockVector{T,F} <: AbstractBlockVector{T,F}
    d::MVector{1,Int64}
    file::F
    LinearBlockVector(file::F,position,::Type{T})  where {F,T} = new{T,F}(MVector{1,Int64}([position]),file)
end
@inline Base.getproperty(cd::LinearBlockVector, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::LinearBlockVector, ::Val{:position_0}) =  :(getfield(cd,:d)[1])
@inline @generated dyncast_get(cd::LinearBlockVector, ::Val{S}) where S = :( getfield(cd, S))
@inline Base.setproperty!(cd::LinearBlockVector, prop::Symbol, val) = dyncast_set(cd,Val(prop),val)
@inline @generated dyncast_set(cd::LinearBlockVector, ::Val{:position_0}, val) = :(getfield(cd,:d)[1]=val)

@inline registerblock(fdfa::LinearBlockVector,::Int64,::Int64) = nothing


# finds the block that contains index, sets file pointer to entry directly after `length` and returns length and current offset of block. 
# If block not yet existent, return 0,0
function find_block(lb::LinearBlockVector{T,F},index::Int) where {T,F}
    position = lb.position_0
    position==0 && (return 0,0)
    seeking = true
    current = 0
    len = 0
    while seeking
        seek(lb.file,position)
        len = read(lb.file,Int64)
        if current+len>=index
            seeking = false
        else
            next = read(lb.file,Int64)
            current += len
            next == 0 && (return 0,0)
            position = next
        end
    end
    return len,current
end

@inline setfirstblock(lb::LinearBlockVector{T,F},index::Int) where {T,F} = (lb.position_0=index)


################################################################################################################

## Vectorized Block-Vectors

################################################################################################################



struct BlockData
    position::Int64
    start::Int64
    length::Int64
end
BlockData() = BlockData(0,0,0)
function Base.write(f::IOStream,b::BlockData)
    write(f,b.position)
    write(f,b.start)
    write(f,b.length)
end
@inline Base.read(f::IOStream,::Type{BlockData}) = BlockData(read(f,Int64),read(f,Int64),read(f,Int64))

fdfa_size(::Type{BlockData}) = 24 # =3*8



struct VectorBlockVector{T,F,DV} <: AbstractBlockVector{T,F}
    d::MVector{1,Int64}
    file::F
    data::DV
    increment::Int64
    VectorBlockVector(file::F,increment,::Type{T})  where {F,T} = new{T,F,Vector{BlockData}}(MVector{1,Int64}([0]),file,Vector{BlockData}(undef,0),increment)
    VectorBlockVector(file::F,increment,::Type{T},dv::DV)  where {F,T,DV} = new{T,F,DV}(MVector{1,Int64}([0]),file,dv,increment)
end
@inline Base.getproperty(cd::VectorBlockVector, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::VectorBlockVector, ::Val{:blockposition}) =  :(getfield(cd,:d)[1])
@inline @generated dyncast_get(cd::VectorBlockVector, ::Val{S}) where S = :( getfield(cd, S))
@inline Base.setproperty!(cd::VectorBlockVector, prop::Symbol, val) = dyncast_set(cd,Val(prop),val)
@inline @generated dyncast_set(cd::VectorBlockVector, ::Val{:blockposition}, val) = :(getfield(cd,:d)[1]=val)

function registerblock(fdfa::VectorBlockVector,position::Int64,size::Int64) 
    fdfa.blockposition =  fdfa.blockposition + 1
    if fdfa.blockposition>length(fdfa.data)
        if typeof(fdfa.data)<:FixedDataFileArray
            resize2(fdfa.data,length(fdfa.data)+fdfa.increment)
        else
            Base.resize!(fdfa.data,length(fdfa.data)+fdfa.increment)
        end
    end
    start = 1
    if fdfa.blockposition>1 
        old_data = fdfa.data[fdfa.blockposition-1]
        start = old_data.start+old_data.length
    end
    fdfa.data[fdfa.blockposition] = BlockData(position,start,size)
end


# finds the block that contains index, sets file pointer to entry directly after `length` and returns length and current offset of block. 
# If block not yet existent, return 0,0
function find_block(hb::VectorBlockVector{T,F},index::Int) where {T,F}
    hb.blockposition==0 && (return 0,0)
    seeking = true
    current = 0
    len = 0
    for i in 1:hb.blockposition
        thisblock = hb.data[i]
        if (thisblock.start+thisblock.length)>index
            seek(hb.file,thisblock.position+fdfa_size(Int64))
            return thisblock.length, thisblock.start-1
        end 
    end
    return 0,0
end

@inline setfirstblock(::VectorBlockVector{T,F},::Int) where {T,F} = nothing



################################################################################################################

## File Arrays

################################################################################################################



# Define the FixedDataFileArray struct for generic type T
struct FixedDataFileArray{T,F,B,BL,V}
    _file::F
    free_allocation::Bool
    position::Int64 # position in file
    buffer::B
    blocklist::BL
    data::MVector{10,Int64}
    visible::V
    #blocks::Int64 # total number of blocks
    #first::Int64 # first block address
    #length::Int64
    #blocktable::Int64 # position of index table for blocks 
    #bufferposition::Int64 # current position within buffer
    #blockposition::Int64 # current buffer_process within block
    #currentblocksize::Int64 # size of the current block
    #offset::Int64
    #blocksize::Int64 # size of any new block (read only)
    #buffersize::Int64 # size of buffers (read only)
    FixedDataFileArray{T}(file, free_allocation,position,_buffer, list,data,visible) where T = new{T,typeof(file),typeof(_buffer),typeof(list),typeof(visible)}(file, free_allocation,position,_buffer, list,data,visible)
end

# array permanently in heap, file IO may vary over time 
const PermanentDataFileArray{T,B,BL,V} = FixedDataFileArray{T,MF,B,BL,V} where {T,B,BL,V,MF<:MutableFile}
# array temporarily in heap with a lifetime shorter than the file IO
const TemporaryDataFileArray{T,B,BL,V} = FixedDataFileArray{T,PF,B,BL,V} where {T,B,BL,PF<:Union{PermanentFile,IO},V}
# array that buffers data in field :buffer
const BufferedDataFileArray{T,F,BL,V} = FixedDataFileArray{T,F,Vector{T},BL,V} where {T,F,BL,V}
# array that reads and writes directly into the file
const NonBufferedDataFileArray{T,F,BL,V} = FixedDataFileArray{T,F,ShortVector{T},BL,V} where {T,F,BL,V}

const NormalFixedDataFileArray{A,B,C,D} = FixedDataFileArray{A,B,C,D,Nothing} where {A,B,C,D}
const TracelessFixedDataFileArray{A,B,C,D} = FixedDataFileArray{A,B,C,D,MVector{1,Int64}} where {A,B,C,D}


@inline Base.length(f::FixedDataFileArray) = f.length

@inline index_in_buffer(fdfa::FDFA,index::Int64) where {T,FDFA<:FixedDataFileArray{T}} = index>fdfa.offset 
function Base.iterate(fdfa::FDFA,start = FileArrayIterator{T}(fdfa)) where {T, FDFA<:FixedDataFileArray{T}}
    increase(start)  
    if start.valid
        return start.current, start
    else
        return nothing
    end
end

function _getindex(fdfa::FDFA,index::Int64) where {T,FDFA<:FixedDataFileArray{T}}
    if index_in_buffer(fdfa,index)
        index-fdfa.offset>fdfa.bufferposition && error("Try to access a field that has not yet been asigned")
        return fdfa.buffer[index-fdfa.offset]
    else
#        block, blockindex = provide_block(fdfa.blocklist,index)
#        return block[blockindex]
        return fdfa.blocklist[index]
    end
end


function _setindex!(fdfa::FDFA,value::T,index::Int) where {T,FDFA<:FixedDataFileArray{T}}
    if fdfa.buffersize==0
        fdfa.blocklist[index] = value
        return value
    end
    if index-fdfa.offset==fdfa.bufferposition+1 # if incrementally last element is assigned
        if fdfa.bufferposition==fdfa.buffersize || (fdfa.currentblocksize!=0 && fdfa.blockposition==fdfa.currentblocksize)
            if fdfa.currentblocksize==0
                    allocate(fdfa,max(fdfa.blocksize,fdfa.buffersize))
                    len,current = find_block(fdfa.blocklist,fdfa.offset+1)
                current!=fdfa.offset && error("not expected")
                fdfa.currentblocksize = len
            end
            moveto(fdfa.blocklist,fdfa.offset+1)
            for i in 1:fdfa.bufferposition
                write(fdfa.file,fdfa.buffer[i])
            end
            fdfa.offset = fdfa.offset + fdfa.bufferposition
            if fdfa.blockposition==fdfa.currentblocksize
                len,current = find_block(fdfa.blocklist,fdfa.offset+1)
                fdfa.currentblocksize = len
                fdfa.blockposition = 0
            end
            fdfa.bufferposition = 0
        end
        fdfa.bufferposition = fdfa.bufferposition + 1
        fdfa.blockposition = fdfa.blockposition + 1
    end
    if index_in_buffer(fdfa,index)
        index-fdfa.offset>fdfa.bufferposition && error("Try to access a field that has not yet been asigned")
        fdfa.buffer[index-fdfa.offset] = value
    else
        fdfa.blocklist[index]=value
    end
    return value
end

@inline function Base.push!(fdfa::FDFA,value::T) where {T,FDFA<:FixedDataFileArray{T}}
    index = fdfa.offset+fdfa.bufferposition+1
    _setindex!(fdfa,value,index)
    return index
end

@inline function Base.push!(fdfa::FDFA,values::AVT) where {T,FDFA<:FixedDataFileArray{T},AVT<:AbstractVector{T}}
    index = 0
    for v in values
        i2 = push!(fdfa,v)
        index = index==0 ? i2 : index 
    end
    return index
end

init(::FDFA) where {FDFA<:NormalFixedDataFileArray} = nothing
exit(::FDFA) where {FDFA<:NormalFixedDataFileArray} = nothing
@inline init(fdfa::FDFA) where {FDFA<:TracelessFixedDataFileArray} = (fdfa.visible[1] = position(fdfa.file))
@inline function exit(fdfa::FDFA) where {FDFA<:TracelessFixedDataFileArray}  
    pos = position(fdfa.file)
    pos!= fdfa.visible[1] && seek(fdfa.file,fdfa.visible[1])
end

@inline Base.getindex(fdfa::FDFA,index::Int) where {T,FDFA<:BufferedDataFileArray{T}} = _getindex(fdfa,index)
@inline Base.getindex(fdfa::FDFA,index::Int) where {T,FDFA<:NonBufferedDataFileArray{T}} = fdfa.blocklist[index]  

@inline Base.setindex!(fdfa::FDFA,value::T,index::Int) where {T,FDFA<:BufferedDataFileArray{T}} = _setindex!(fdfa,value,index)
@inline function Base.setindex!(fdfa::FDFA,value::T,index::Int) where {T,FDFA<:NonBufferedDataFileArray{T}}
    init(fdfa)
    b = index_in_buffer(fdfa,index)
        _setindex!(fdfa,value,index)
        b && _setindex!(fdfa,value,index+1)
        exit(fdfa)
end



@inline function write_data_position(cd::TemporaryDataFileArray,i,val)
    getfield(cd, :data)[i] = val
    init(cd)
    seek(cd.file,cd.position + (i-1)*sizeof(Int64))
    write(cd.file,getfield(cd, :data)[i])
    exit(cd)
end
@inline write_data_position(cd::PermanentDataFileArray,i,val) = getfield(cd, :data)[i] = val

@inline Base.getproperty(cd::FixedDataFileArray, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::FixedDataFileArray, ::Val{:file}) =  :(getfield(cd,:_file))
@inline @generated dyncast_get(cd::FixedDataFileArray, ::Val{:blocks}) =  :(getfield(cd,:data)[1])
@inline @generated dyncast_get(cd::FixedDataFileArray, ::Val{:first}) =  :(getfield(cd,:data)[2])
@inline @generated dyncast_get(cd::FixedDataFileArray, ::Val{:length}) =  :(getfield(cd,:data)[3])
@inline @generated dyncast_get(cd::FixedDataFileArray, ::Val{:blocktable}) =  :(getfield(cd,:data)[4])
@inline @generated dyncast_get(cd::FixedDataFileArray, ::Val{:bufferposition}) =  :(getfield(cd,:data)[5])
@inline @generated dyncast_get(cd::FixedDataFileArray, ::Val{:blockposition}) =  :(getfield(cd,:data)[6])
@inline @generated dyncast_get(cd::FixedDataFileArray, ::Val{:currentblocksize}) =  :(getfield(cd,:data)[7])
@inline @generated dyncast_get(cd::FixedDataFileArray, ::Val{:offset}) =  :(getfield(cd,:data)[8])
@inline @generated dyncast_get(cd::FixedDataFileArray, ::Val{:blocksize}) =  :(getfield(cd,:data)[9])
@inline @generated dyncast_get(cd::FixedDataFileArray, ::Val{:buffersize}) =  :(getfield(cd,:data)[10])
    #blocktable::Int64 # position of index table for blocks 
    #bufferposition::Int64 # current position within buffer
    #blockposition::Int64 # current buffer_process within block
    #currentblocksize::Int64 # size of the current block
    #offset::Int64
    #blocksize::Int64 # size of any new block (read only)
    #buffersize::Int64 # size of buffers (read only)

@inline @generated dyncast_get(cd::FixedDataFileArray, ::Val{S}) where S = :( getfield(cd, S))
@inline Base.setproperty!(cd::FixedDataFileArray, prop::Symbol, val) = dyncast_set(cd,Val(prop),val)
@inline @generated dyncast_set(cd::FixedDataFileArray, ::Val{:blocks}, val) = :(write_data_position(cd,1,val))
@inline @generated dyncast_set(cd::FixedDataFileArray, ::Val{:first}, val) = :(write_data_position(cd,2,val))
@inline @generated dyncast_set(cd::FixedDataFileArray, ::Val{:length}, val) = :(write_data_position(cd,3,val))
@inline @generated dyncast_set(cd::FixedDataFileArray, ::Val{:blocktable}, val) = :(write_data_position(cd,4,val))
@inline @generated dyncast_set(cd::FixedDataFileArray, ::Val{:bufferposition}, val) = :(write_data_position(cd,5,val))
@inline @generated dyncast_set(cd::FixedDataFileArray, ::Val{:blockposition}, val) = :(write_data_position(cd,6,val))
@inline @generated dyncast_set(cd::FixedDataFileArray, ::Val{:currentblocksize}, val) = :(write_data_position(cd,7,val))
@inline @generated dyncast_set(cd::FixedDataFileArray, ::Val{:offset}, val) = :(write_data_position(cd,8,val))
@inline @generated dyncast_set(cd::FixedDataFileArray, ::Val{:blocksize}, val) = :(write_data_position(cd,9,val))
@inline @generated dyncast_set(cd::FixedDataFileArray, ::Val{:buffersize}, val) = :(write_data_position(cd,10,val))


struct Slow
    Slow() = new()
    Slow(::Int64) = new()
end
@inline Base.getproperty(cd::Slow, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(::Slow, ::Val{:blocksize}) =  :(return 0)
@inline @generated dyncast_get(cd::Slow, ::Val{S}) where S = :( getfield(cd, S))


struct Fast{T<:StaticBool}
    blocksize::Int64
    buffersize::Int64
    buffered::T
end

struct HeapFast
    blocksize::Int64
end

# Constructor in case array already exists on harddrive or needs to be created
function FixedDataFileArray{T}(file::AF, position::Int64,free_allocation::Bool,buffered::StaticBool,blocklist,_buffersize::Int64,blocksize::Int64,v=nothing) where {AF<:Union{AbstractFile,IO},T}
    _buffersize = max(1,_buffersize)
    buffer(::StaticTrue) = _buffersize, Vector{T}(undef,_buffersize)
    buffer(::StaticFalse) = 1, ShortVector{T}()
    buffersize, _buffer = buffer( StaticBool(buffered))
    adjust_block(buffer::Int, block::Int) = block % buffer == 0 ? block : ((block ÷ buffer) + 1) * buffer
    blocksize = adjust_block(buffersize,blocksize)
    if position==0
        position = allocate(file,sizeof(fieldtype(FixedDataFileArray,:data)),free_allocation)
        ldata = size(fieldtype(FixedDataFileArray,:data))[1]
        write_Int64_position(file,position,ldata-1,[blocksize,buffersize])
    end
    # Go to the specified position
    seek(file, position)
    
    # Read the blocks and first value from the file
    data = read(file,fieldtype(FixedDataFileArray,:data))
     
    # Return the new instance
    #return FixedDataFileArray(file, free_allocation,position,_buffer, list(_blocklist,_heapblocklist,data[5]),data)
    getblocklist(h::HeapFast) = false, VectorBlockVector(file,h.blocksize,T)
    getblocklist(::Slow) = false, LinearBlockVector(file,data[2],T)
    function getblocklist(f::Fast)
        # recall data[4] == blocktable
        my_data = FixedDataFileArray{BlockData}(file,data[4],true,f.buffered,Slow(),f.buffersize,f.blocksize,MVector{1,Int64}([0]))
        data[4] = my_data.position
        return data[4] != my_data.position, VectorBlockVector(file,f.blocksize,T,my_data)
    end
    changeblocktable, myblocklist = getblocklist(blocklist)
    fdfa = FixedDataFileArray{T}(file, free_allocation,position,_buffer, myblocklist,data,v)
    changeblocktable && (fdfa.blocktable = myblocklist.position)
    return fdfa
end

function allocate(fdfa::FixedDataFileArray{T},size::Int64) where T
    block_size = size*fdfa_size(T)  
    position = allocate(fdfa.file,block_size,fdfa.free_allocation,data=[size,0])
    if fdfa.length>0 
        find_block(fdfa.blocklist,fdfa.length)
        write(fdfa.file,position)
    else
        fdfa.first = position
        setfirstblock( fdfa.blocklist,position)
        fdfa.currentblocksize = size
    end
    fdfa.length = fdfa.length + size
    fdfa.blocks = fdfa.blocks + 1
    registerblock(fdfa.blocklist,position,size)
end
@inline Base.resize!(fdfa::FDFA,size::Int64) where {FDFA<:FixedDataFileArray} = allocate(fdfa,size-fdfa.length)
@inline resize2(fdfa::FDFA,size::Int64) where {FDFA<:FixedDataFileArray} = allocate(fdfa,size-fdfa.length)


#############################################################################################################################

## 64 Bit Container

#############################################################################################################################

struct Container64{F,B,BL,V}
    data::FixedDataFileArray{UInt64,F,B,BL,V}
    inds::MVector{3,Int64}
    UINTBuffer::Vector{UInt64}
    INTBuffer::Vector{Int64}
    FloatBuffer::Vector{Float64}
    Container64(d::FixedDataFileArray{UInt64,F,B,BL,V}) where {F,B,BL,V} = new{F,B,BL,V}(d,MVector{3,Int64}([0,0,0]),Vector{UInt64}(undef,0),Vector{Int64}(undef,0),Vector{Float64}(undef,0))
    function Container64(file::AF, position::Int64,free_allocation::Bool,buffered::StaticBool,blocklist,_buffersize::Int64,blocksize::Int64,v=nothing) where {AF<:Union{AbstractFile,IO}}
        Container64(FixedDataFileArray{UInt64}(file, position,free_allocation,buffered,blocklist,_buffersize,blocksize,v) )
    end
end

@inline Base.getproperty(cd::Container64, prop::Symbol) = dyncast_get(cd,Val(prop))
@inline @generated dyncast_get(cd::Container64, ::Val{:index}) =  :(getfield(cd,:inds)[1])
@inline @generated dyncast_get(cd::Container64, ::Val{:position}) =  :(getfield(cd,:inds)[2])
@inline @generated dyncast_get(cd::Container64, ::Val{:remaining_block_length}) =  :(getfield(cd,:inds)[3])
@inline @generated dyncast_get(cd::Container64, ::Val{:UINTBuffer}) =  :(getfield(cd,:UINTBuffer))
@inline @generated dyncast_get(cd::Container64, ::Val{:INTBuffer}) =  :(getfield(cd,:INTBuffer))
@inline @generated dyncast_get(cd::Container64, ::Val{:FloatBuffer}) =  :(getfield(cd,:FloatBuffer))
@inline @generated dyncast_get(cd::Container64, ::Val{S}) where S =  :(getproperty(getfield(cd,:data),S))

@inline Base.setproperty!(cd::Container64, prop::Symbol, val) = dyncast_set(cd,Val(prop),val)
@inline @generated dyncast_set(cd::Container64, ::Val{:index}, val) = :(getfield(cd,:inds)[1]=val)
@inline @generated dyncast_set(cd::Container64, ::Val{:position}, val) = :(getfield(cd,:inds)[2]=val)
@inline @generated dyncast_set(cd::Container64, ::Val{:remaining_block_length}, val) = :(getfield(cd,:inds)[3]=val)
@inline @generated dyncast_set(cd::Container64, ::Val{S}, val) where S = :(setproperty!(getfield(cd,:data),S,val))

@inline Base.getindex(c::Container64{F, B, BL, V}, index) where {F, B, BL, V} = getindex(getfield(c, :data), index)
@inline Base.setindex!(c::Container64{F, B, BL, V}, value, index) where {F, B, BL, V} = setindex!(getfield(c, :data), value, index)
@inline Base.length(c::Container64{F, B, BL, V}) where {F, B, BL, V} = length(getfield(c, :data))
@inline Base.iterate(c::Container64{F, B, BL, V}) where {F, B, BL, V} = iterate(getfield(c, :data))
@inline Base.iterate(c::Container64{F, B, BL, V}, state) where {F, B, BL, V} = iterate(getfield(c, :data), state)


const BIT64 = Union{Int64,Float64,UInt64}

@inline Base.push!(c::C64,data::BIT64) where {C64<:Container64} = push!(getfield(c, :data),reinterpret(UInt64,data))

@inline function Base.push!(c::C64,data::AVB) where {C64<:Container64,B<:BIT64,AVB<:AbstractVector{B}}
    len = length(data)
    if length(c.UINTBuffer)<len 
        resize!(c.UINTBuffer,len)
    end
    for i in 1:len
        c.UINTBuffer[i] = reinterpret(UInt64,data[i])
    end
    return push!(getfield(c, :data),view(c.UINTBuffer,1:len))
end

@inline function Base.push!(c::C64,data::P) where {C64<:Container64,AV1<:AbstractVector{Int64},AV2<:AbstractVector{Float64},P<:Pair{AV1,AV2}}
    pos = push!(c,length(data[1]))
    push!(c,data[1])
    push!(c,length(data[2]))
    push!(c,data[2])
    return pos
end

function moveto(c::C64,index::Int) where {C64<:Container64}
    if index>c.offset 
        c.index = index
        return
    end
    pos = position(c.file)
    delta = 0
    if pos!=c.position || index>=(c.index+c.remaining_block_length) || index<c.index
        rbl, current = find_block(c.blocklist,index)
        c.remaining_block_length = rbl
        println("here $(c.remaining_block_length), $(find_block(c.blocklist,index))")
        read(c.file,UInt64)
        delta = index-current-1
    else
        delta = index-c.index
    end
    for i in 1:delta
        read(c.file,UInt64)
    end
    c.index = index
    c.remaining_block_length = c.remaining_block_length - delta
    c.position = position(c.file)
end

@inline function Base.read(c::C64,::Type{T}) where {C64<:Container64,T<:BIT64} 
    if c.index<=c.offset
        if c.remaining_block_length<=0
            moveto(c,c.index)
        end
        c.index = c.index + 1
        c.remaining_block_length = c.remaining_block_length - 1
        c.position = c.position + 8
        return  reinterpret(T,read(c.file,UInt64))
    else
        c.index = c.index + 1
        return reinterpret(T,getfield(c,:data)[c.index-1])
    end
end

@inline function Base.read(c::C64,::Type{T},index::Int64) where {C64<:Container64,T} 
    moveto(c,index)
    read(c,T)
end

function readblock(c::C64,len::Int64,index::Int64,::Type{T}) where {C64<:Container64,T} 
    moveto(c,index)
    readblock(c,len,T)
end
function readblock(c::C64,len::Int64) where {C64<:Container64} 
    len>length(c.UINTBuffer) && resize!(c.UINTBuffer,len)
    for i in 1:len
        c.UINTBuffer[i] = read(c,UInt64)
    end
    return view(c.UINTBuffer,1:len)
end

function readblock(c::C64,len::Int64,::Type{Int64}) where {C64<:Container64} 
    len>length(c.INTBuffer) && resize!(c.INTBuffer,len)
    for i in 1:len
        c.INTBuffer[i] = reinterpret(Int64,read(c,UInt64))
    end
    return view(c.INTBuffer,1:len)
end

function readblock(c::C64,len::Int64,::Type{Float64}) where {C64<:Container64} 
    len>length(c.FloatBuffer) && resize!(c.FloatBuffer,len)
    for i in 1:len
        c.FloatBuffer[i] = reinterpret(Float64,read(c,UInt64))
    end
    return view(c.FloatBuffer,1:len)
end
    

struct HVVertex
end

@inline function Base.read(c::C64,::Type{HVVertex}) where {C64<:Container64} 
    len = read(c,Int64)
    b1 = readblock(c,len,Int64)
    len = read(c,Int64)
    b2 = readblock(c,len,Float64)
    return b1, b2
end



###############################################################################################################################

## HVFile & Container

###############################################################################################################################

const internalfilenumber = MVector{1, Int64}([0])

function current_filename()
    c = internalfilenumber[1]
    internalfilenumber[1] += 1
    return "HV" * lpad(c, 6, '0') * ".hvd"
end




mutable struct HVFile{C64<:Container64} <: AbstractFile
    filename::String
    file::MutableFile{IOStream}
    open::Bool
    container::C64
    function HVFile(filename::String, _buffersize::Int64=2^17, blocksize::Int64=2^20;overwrite=false)
        if !isfile(filename) || overwrite
            f = open(filename, "w")
            write(f, Int64(1))  # Schreibe eine initialisierende Int64
            close(f)
            f = open(filename, "a+")
            seekend(f)
            println(position(f))
            c = Container64(MutableFile(f), 0, true, StaticTrue(), HeapFast(5), _buffersize, blocksize)
            seek(f,0)
            for i in 1:11
                print("$(read(f,Int64)) - ")
            end
            println()
            println(c.position)
            close(f)
        end
        f = open(filename, "a+")
        mf = MutableFile(f)
        c = Container64(mf, 8, true, StaticTrue(), HeapFast(5), _buffersize, blocksize)
        ret = new{typeof(c)}(filename, mf, false, c)
        close(f)
        return ret
    end
    HVFile(_buffersize::Int64=2^17, blocksize::Int64=2^20) = HVFile(current_filename(),_buffersize,blocksize)
end

function Base.open(func::Function,f::HVF) where {HVF<:HVFile}
    was_open = f.open
    try
        if !f.open
            f.file.file = open(f.filename, "a+")
            f.open = true
        end
        func(f)
    catch err
        if !was_open
            f.open = false
            close(f.file)
        end
        rethrow(err)
    finally
        if !was_open
            f.open = false
            close(f.file)
        end
    end
end


# Just a buffer if we have no database....

struct NoFile end
@inline Base.open(func::Function,vg::NoFile) = func(vg)


##### sample usage copied
#=
# Öffne die Datei im Schreibmodus und schreibe Int64(1)
open("test.dat", "w") do io
    write(io, Int64(1))
end

# Öffne die Datei erneut mit einem `do` Block für weitere Operationen
open("test.dat", "a+") do io
    array = FixedDataFileArray{Int32}(io,0,true,statictrue,Fast(10,5,statictrue),5,10)
    #array = FixedDataFileArray{Int64}(io,0,true,statictrue,Slow(),5,10)
    println(array.position)
        # Liest Int64-Werte und gibt sie getrennt durch Kommas aus
    for i in 1:200
        array[i] = Int32(i)
    end 
#end
#open("test.dat", "r") do io
seek(io,0)
first = true
    while !eof(io)
        value = read(io, Int64)
        if first
            print(value)
            first = false
        else
            print(",", value)
        end
    end
    println("")
    for i in 1:200
        print("$(array[i]), ")
    end
# Hier können weitere Operationen eingefügt werden
    # Zum Beispiel: Lesen und verarbeiten der Daten
end

=#
