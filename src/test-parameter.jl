#using Iterators
function para(;kwargs...)
    inds=collect(1:(length(kwargs)))
    print("first : ")
    println(first(zip(keys(kwargs),inds)))
    println(keys(kwargs),values(kwargs))
    println("erster: $(kwargs[2])")
end

tupi=(alpha=1, gamma="Hallo")

println(keys(tupi))

for k in tupi
    println("$k")
end

para(;tupi...)

function f_vcat(;kwargs...)
    return x->vcat()
end

vcat(Float64[],[1,2,3],1)

map(x->x^2,(alpha=1,beta=3.1,eta=2))

#typeof((alpha=1,beta=2))==typeof((alpha=3,beta=20))

vcat(Float64[],(alpha=1,beta=3.1,eta=2)...)

using SpecialFunctions

g(x)=map(f->f(x),(alpha=exp,beta=sin,eta=y->y^2))

g(3)
length(1.5)

function emptyf(;kw...)
    println(kw...)
    emptyf2(hallo=57;values(kw)...)
end


function emptyf2(;kw...)
    println("here",values(kw),keys(kw))
    map(x->x^2+1,values(kw))
end

emptyf(alpha=1,beta=2)#;tupi...)
    
#println(Base.Pairs(tupi,keys(tupi)))

tuppiii=(1,"a",3.2)

println(tuppiii[3])
println(isdefined(tuppiii,4))

function unbox((a,b,c))
    println(c)
end
unbox(tuppiii)
println(typeof(tuppiii)<:Tuple)
println(typeof([1,2])<:Tuple)
println(typeof(tupi)<:NamedTuple)
println(typeof(tuppiii)<:NamedTuple)
println(typeof(zeros(Int,2))<:Vector{Int})

a=[[1.0,2.3],[1.3]]
b=[[1.3,4.3,5.6]]

append!(b,a)

a[1][2]=6.7

println(b)

f=(;x...)->(y=x[:alpha]^2;
    return y+x[:beta])

println(f(alpha=3,beta=8,gamma=9))

println(typeof(:alpha))

data=rand(2,4)
data[:,2]=[1.0,2.0]

println(data[:,2])

ddd=Dict(1=>3.4,7=>"Hallo")

typeof(ddd)

using InteractiveUtils

foo(a_type::Type, an_interface::Symbol) = an_interface ∈ [i.name for i in methodswith(a_type,supertypes=true)]

#println(foo(typeof(ddd),:getindex))
#println(foo(typeof(emptyf),:getindex))

ff=[1,2,3,emptyf]

println(typeof(ff))
k=ff[4]
println(typeof(k))

iii= [3,4,7]
zz=zeros(Float64,10)

#println(map(k->(zz[k]=1.0),iii))

broadcast(k->(zz[k]=1.0),iii)

println(zz)

println(typeof([1,2])<:Function)

a=(4,"h")
b=(5.6,3)
println((a...,b...))

function _emptyf(x;kw...)
    println(x,values(kw))
end
_emptyf("teest ",alpha=10,gamma="s")
ffff2(;kw...)=_emptyf("100 $(kw[:alpha]) : ";kw...)

ffff2(alpha=1,beta=10.4)

