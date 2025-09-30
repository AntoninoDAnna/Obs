@doc raw"""
     sym_der(v,bnd::Boundary)

Compute the symmetryc derivative of `v` according to the boundary condition in `bnd`
See also [`Boundary`](@ref)
"""
function sym_der(v,bnd::Boundary)
    res = similar(v)
    res[2:end-1] = 0.5 .* (v[3:end].-v[1:end-2])
    if bnd == open
        res[1] = v[2] - v[1]
        res[end] = v[end] - v[end-1]
    elseif bnd == periodic
        res[1] = 0.5*(v[2]-v[end])
        res[2] = 0.5*(v[1]-v[end-1])
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
function sym_source(corr, y0, parity, bnd::Boundary)
    T = lastindex(corr)
    if bnd == open
        ie,is =  (T>2y0) ? (2y0-1,1) : (T,2y0-T)
        res = 0.5.*(corr[y0:ie] .+ parity*corr[y0:-1:is])
    elseif bnd ==periodic
        res = zeros(typeof(corr[1]),div(T,2))
        for i in eachindex(res)
            res[i] = 0.5*( corr[(y0+i-1)%T] + corr[(y0-i+1)%T] )
        end
    end
    return res
end

function sym_source(corr::juobs.Corr,y0,parity,bnd::Boundary)
    obs = if bnd == open
       sym_source(corr.obs[2:end-1],y0,parity,bnd)
    elseif bnd == periodic
        sym_source(corr.obs,y0,parity,bnd)
    end
    return juobs.Corr(obs,corr.kappa,corr.mu,corr.gamma,0,corr.theta1,corr.theta2)
end

sym_source(corr::juobs.Corr, parity, bnd::Boundary) = sym_source(corr,corr.y0,parity,bnd)


@doc """
     check_corr(c::juobs.Corr...; flag::Check_flag)

check that all correlator are compatible. Use flag to ignore specific fields

See also [`Check_flag`](@ref)
"""
function check_corr(c::juobs.Corr...; flag::Check_flag)
    fields =let
        aux =(flag .& instances(Check_flag)[2:end] .==no_flag) |> collect # remove no_flag and no_thetas
        [:gamma, :obs, :kappa, :mu, :y0, :theta1, :theta2][aux]
    end
    if :obs in fields
        x = getfield.(c,:obs)
        n = length(x[1])
        if any(length.(x[2:end]).!=n)
            error("[check_corr] No compatible correlators. Fields :obs have different lengths")
        end
        filter!(x->x ≠ :obs, fields)
    end
    for f in fields
        x = getfield.(c,f)
        if any(x .!= [x[1]])
            error("[check_corr] No compatible correlators. Fields $f are different")
        end
    end
end

@doc """
     average_corr(x::juobs.Corr...; flag::Check_flag = no_gamma)

It average the correlator `x...`. It first check for compatibility among the
correlator according to `flag`. It is meant to average correlator with vector
and tensor currents and `flag = no_gamma`.

## Output Specifics

The function returns a `juobs.Corr` with gamma structure made by appending
the gamma structure of the original correlator, so that if we average
three correlator with gamma structure `G1`, `G2` and `G3`, the returned
correlator has gamma structure `G1,G2,G3`.

#TODO: if `flags` is something different that `no_gamma`, the correlators will
be averaged correctly but other informations will be lost (`:kappa`, `:mu`, ecc...)

See also [`check_corr`](@ref), [`Check_flag`](@ref)
"""
function average_corr(x::juobs.Corr...;flag::Check_flag = no_gamma)
    check_corr(x..., flag = flag)
    Nc = length(x)
    obs = getfield.(x,:obs)
    gamma = getfield.(x,:gamma)
    G1 = if all(x->x[1]==gamma[1][1], gamma[2:end])
        gamma[1][1]
    else
        join([g[1] for g in gamma],",")
    end
    G2 = if all(x->x[2]==gamma[1][2], gamma[2:end])
        gamma[1][2]
    else
        join([g[2] for g in gamma],",")
    end

    mean = reduce(+,obs)/Nc
    return  juobs.Corr(mean,x[1].kappa,x[1].mu, [G1,G2],x[1].y0,x[1].theta1,x[2].theta2)
end
