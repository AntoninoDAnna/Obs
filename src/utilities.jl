@doc raw"""
     sym_der(v,bnd::Boundary)

Compute the symmetryc derivative of `v` according to the boundary condition in `bnd`
See also [`Boundary`](@ref)
"""
function sym_der(v,::OBC)
    res = similar(v)
    for i in 2:lastindex(res)-1
        res[i] = 0.5 .* (v[i+1].-v[i-1])
    end
    res[1] = v[2] - v[1]
    res[end] = v[end] - v[end-1]
    return res
end

function sym_der(v,::PBC)
    res = similar(v)
    T = length(res)
    for i in eachindex(res)
        res[i] = 0.5 .* (v[i%T+1].-v[(i-T-2)%T+1])
    end
    return res
end


@doc raw"""
    sym_source(corr, y0, parity, bnd::Boundary = open)
    sym_source(corr::juobs.Corr [,y0], parity, bnd::Boundary = open)

symmetrize the correlator with respect to the source position `y0`

## Arguments
  - `y0` source position in Julia indexing.
  - `parity` parity of the correlator. should be `\pm 1`
  - `bnd` the boundary condition of the correlator

## Return types
  - (1) returns a vector of the same type as `corr`
  - (2) returns a `juobs.Corr`

See also [`Boundary`](@ref)
"""
function sym_source(corr, y0, parity, ::OBC)
    T = lastindex(corr)
    is, ie = if T == 2y0
        1,T-1
    elseif T<2y0
        2y0-T, T
    else
        1,2y0-1
    end
    res = 0.5.*(corr[y0:ie] .+ parity*corr[y0:-1:is])
    return res
end

function sym_source(corr,y0,parity,::PBC)
    res = zeros(typeof(corr[1]),div(T,2))
    for i in eachindex(res)
        res[i] = 0.5*( corr[(y0+i-1)%T] + corr[(y0-i+1)%T] )
    end
    return res
end

sym_source(corr,y0,parity) = sym_source(corr,y0,parity,OBC())

sym_source(corr::AbstractCorr,y0,parity,::OBC) = sym_source(corr.obs[2:end-1],y0,parity)
sym_source(corr::AbstractCorr,y0,parity,::PBC) = sym_source(corr.obs,y0,parity)

sym_source(corr::AbstractCorr, parity,bnd::Union{OBC,PBC}) =
    sym_source(corr,ObsIO.src(corr),parity,bnd)

sym_source(corr::AbstractCorr,parity) = sym_source(corr,ObsIO.src(corr),parity,OBC())


function get_average_point(C,n)
    P = C[1].points[n]
    any(P.gamma != c.points[n].gamma for c in C ) || return P
    g = ntuple(i->C[i].points[n].gamma,length(C))
    if all(x in g for x in (ObsIO.G1,ObsIO.G2,ObsIO.G3))
        return ObsIO.__update__(P,gamma = ObsIO.Gi)
    end
    if all(x in g for x in (ObsIO.G0G1,ObsIO.G0G2,ObsIO.G0G3))
        return ObsIO.__update__(P,gamma = ObsIO.G0Gi)
    end
    if all(x in g for x in (ObsIO.G1G5,ObsIO.G2G5, ObsIO.G3G5))
        return ObsIO.__update__(P,gamma = ObsIO.GiG5)
    end
    if all(x in g for x in (ObsIO.G1G2,ObsIO.G2G3, ObsIO.G1G3))
        return ObsIO.__update__(P,gamma = ObsIO.GiGj)
    end
    error("cannot average these points")
end

function get_averaged_prop(C,pts,n)
    pr = C[1].propagators[n]
    for c in C, f in  (:k,:mu,:pF,:seq_prop,:theta)
        if !isequal(getfield(pr,f),getfield(c.propagators[n],f))
            error("cannot average these propagators")
        end
    end
    if any(pr.src != c.propagators[n].src for c in C)
        _p = findfirst(isequal(pr.src.x0 ,p.x0) for p in pts)
        return ObsIO.__update__(pr,src = pts[_p])
    end
    if any(pr.snk != c.propagators[n].snk for c in C)
        _p = findfirst(isequal(pr.snk.x0, p.x0) for p in pts)
        return ObsIO.__update__(pr,snk = pts[_p])
    end
    return pr
end

@doc """
     average_corr(x::ObsIO.Corr{N}) where N

It average the correlator `x...`.

## Output Specifics

The function returns a `Corr{N}` with gamma structure replace as follow:
  - all gammas are equal -> keep the current structure
  - G1,G2,G3 -> Gi
  - G0G1,G0G2,G0G3 -> G0Gi
  - G5G1,G5G2,G5G3 -> G5Gi
  - G1G2, G1G3, G2G3 -> GiGj
"""
function average_corr(C::ObsIO.Corr{N}...) where N
    Nc = length(C)
    f(n)  = get_average_point(C,n)
    pts = ntuple(f,N)
    g(n) = get_averaged_prop(C,pts,n)
    prp = ntuple(g,N)
    obs = reduce(+,getfield.(C,:obs))./Nc
    return ObsIO.Corr(obs,pts,prp)
end


@doc raw"""
    plat_av(obs::T where {T<:AbstractVector} ,W::AbstractVecOrMat{<:Real})

Compute the plateau average of `obs`. If `W` is an `AbstractVector`, then the
plateau average correspond to the analytic solution of an uncorrelated
fit with a constant function, if `W` is a `AbstractMatrix`, then the plateau
average correspond to the analytic solution of a correlated fit with a
constant function.

"""
plat_av(obs::T where {T<:AbstractVector}, W::AbstractMatrix{<:Real}) =
    sum(obs[i]*W[i,j] for i in axes(W,1), j in axes(W,2))/sum(W)

plat_av(obs::T where {T<:AbstractVector}, W::AbstractVector{<:Real}) =
    sum(W.*obs)/sum(W)
