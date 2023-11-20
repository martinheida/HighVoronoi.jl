########################################################################################################

# calculate a determinant iteratively according to the Leibnitz Minor method.

########################################################################################################

struct Minors
    HASH::Matrix{Int64}
    data::Vector{Vector{Float64}}
    Index::Vector{Int64}
    dimension::Int64
    _buffer::Vector{Vector{Int64}}
    function Minors(dim)
        H=hash_matrix(dim)
        d=Vector{Vector{Float64}}(undef,dim)
        b=Vector{Vector{Int64}}(undef,dim)
        for col in 1:dim 
            d[col]=zeros(Float64,hash_column(H,col,dim)) 
            b[col]=zeros(Int64,col)
        end
        Ind=zeros(Int64,dim)
        new(H,d,Ind,dim,b)
    end
end

function hash_matrix(dim)
    HASH=zeros(Int64,dim,dim)
    for j in 1:dim HASH[1,j]=1 end
    for i in 2:dim
        for j in 1:(dim-(i-1))
            for l in (j+1):(dim-(i-1)+1)
                HASH[i,j]+=HASH[i-1,l]
            end
        end
    end
    return HASH
end

function vor_minor_determinant(minors,vectors)
    for i in 1:(minors.dimension)
        k_minor(minors,i,vectors[i])
        #println(all_determinants[i])
    end
    
    return (minors.data[end])[1]
end



function hash_column(hash,col,dim)
    _sum=0
    for i in 1:dim _sum+=hash[col,i] end
    return _sum
end

function k_minor(minors::Minors,k,V)
    #offset=0
    dim=minors.dimension
    #println(k)
    if k==1
        data=minors.data[1]
        for i in 1:dim data[i]=V[i] end
    else
        #for i in 1:k offset+=hash_column(minors.HASH,i,dim) end
        _minor=minors._buffer[k] # This represents a minor of k rows (sorted) and the last k columns
        for j in 1:k _minor[j]=j end # minor has the numbers (1,2,3,....,k)
        index=1
        #println(_minor)
        while true
            result=0.0
            #println(_minor)
            DATAk_m_one=(minors.data[k-1])
            sign=(-1)
            for j in 1:k 
                sign*=(-1) # same as factor (-1)^(j+1)
                result+=V[_minor[j]]*(DATAk_m_one[get_minor(minors,_minor,j,k-1)])*sign 
                #println("V($(_minor[j]))*M($_minor,$j)=$(V[_minor[j]])*$(((minors.data[k-1])[get_minor(minors,_minor,j,k-1)]))")
            end
            #(minors.data[k])[get_minor(minors,_minor,0,k)]=result
            (minors.data[k])[index]=result # same result, but faster....
            index+=1
            if !next_minor(_minor,dim) break end
        end
    end
end

function next_minor(minor,dim)
    k=length(minor)
    column=k
    while true
        if column<=0 return false end
        minor[column]+=1
        if minor[column]<=dim-k+column
            for _cc in (column+1):k
                minor[_cc]=minor[_cc-1]+1
            end
            return true
        else
            column-=1
        end
    end
end

function get_minor(minors::Minors,indeces,skip,top_col)
    # minors: list of minors, indeces: array storing current indeces, 
    # _cancel: position in ideces which is not considered, _length: length of sublock, i.e. length(indeces)-1
    index=1
    storage=0
    col=top_col
    old_row=0
    while col>0
        if index==skip index+=1 end
        for i in (old_row+1):indeces[index]-1 
            storage+=minors.HASH[col,i]
        end
        old_row=indeces[index]
        storage+= col==1 ? 1 : 0
        col-=1
        index+=1
    end
    #println("ind: $indeces, skip: $skip, storage: $storage")
    return storage
end

