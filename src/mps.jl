mutable struct MPS
    sites::Array{Site, 1}
    infinite::Bool
end

MPS(sites; infinite=true) = MPS(sites, infinite)

Base.length(mps::MPS) = length(mps.sites)
Base.getindex(mps::MPS, i::Int) = mps.sites[i]
# TODO is it desired to inherit the 'infinite' status?
Base.getindex(mps::MPS, r::UnitRange) = MPS(mps.sites[r], mps.infinite)
Base.lastindex(mps::MPS) = length(mps.sites)

function Base.setindex!(mps1::MPS, mps2::MPS, r::UnitRange)
    for (i, ri) in enumerate(r)
        mps1.sites[ri] = mps2.sites[i]
    end
end

function Base.iterate(mps::MPS)
    if length(mps) > 0
        return (mps.sites[1], 1)
    else
        return nothing
    end
end

function Base.iterate(mps::MPS, i)
    if i < length(mps.sites)
        return (mps.sites[i+1], i+1)
    else
        return nothing
    end
end

isinfinite(mps::MPS) = mps.infinite

function get_Θ(mps::MPS)
    Θ = if isinfinite(mps)
        get_Λ(mps[end])
    else
        Diagonal([1.0])
    end
    χL, χR = size(Θ)
    Θ = collect(reshape(Θ, (χL, 1, χR)))
    for (i, site) in enumerate(mps)
        @tensor begin
            Θ[α, i, j, β] := Θ[α, i, γ] * get_A(site)[γ, j, β]
        end
        χL, di, dj, χR = size(Θ)
        Θ = reshape(Θ, (χL, di*dj, χR))
    end
    return Θ
end

function get_Θ(mps::MPS, r::UnitRange)
    a, b = r.start, r.stop
    Θ = get_left_Λ(mps, a)
    χL, χR = size(Θ)
    Θ = collect(reshape(Θ, (χL, 1, χR)))
    for (i, site) in enumerate(mps[a:b])
        @tensor begin
            Θ[α, i, j, β] := Θ[α, i, γ] * get_A(site)[γ, j, β]
        end
        χL, di, dj, χR = size(Θ)
        Θ = reshape(Θ, (χL, di*dj, χR))
    end
    return Θ
end

function distance(mps1::MPS, mps2::MPS)
    return norm(vec(get_Θ(mps1)) .- vec(get_Θ(mps2)))
end

function get_bond_dims(mps::MPS)
    return [get_right_bond_dim(site) for site in mps]
end

function get_entanglements(mps::MPS)
    return [get_entanglement(site) for site in mps]
end

function LinearAlgebra.norm(mps::MPS)
    S = get_Θ(mps)
    @tensor begin
        n[:] := S[α, i, β] * conj(S)[α, i, β]
    end
    return real(n)
end

function expectation(mps::MPS, O)
    S = get_Θ(mps)
    @tensor begin
        e[:] := S[α, i, β] * O[i, j] * conj(S)[α, j, β]
    end
    return real(e) / norm(mps)
end

function get_left_Λ(mps::MPS, i)
    return if i == 1
        if isinfinite(mps)
            get_Λ(mps[end])
        else
            Diagonal([1.0])
        end
    else
        get_Λ(mps[i-1])
    end
end

function shuffle!(mps::MPS)
    n = length(mps)
    half = div(n, 2)
    mps.sites = vcat(mps.sites[(half+1):end], mps.sites[1:half])
    return mps
end

function decompose(Θ, χ_max, trim)
    # TODO what is a good condition here?
    F = if (size(Θ, 1) < 20) && (size(Θ, 2) < 20)
        svd(Θ)
    else
        # @show size(Θ), min(χ_max, size(Θ, 1)-2)
        svds(Θ, nsv=min(χ_max, size(Θ, 1)-2, size(Θ, 2)-2))[1]
    end
    u, s, v = F.U, F.S, F.Vt
    gt = count(λ -> λ >= trim, s)
    χ = min(length(s), χ_max, gt)
    u = u[:, 1:χ]
    s = s[1:χ]
    n = norm(s)
    s /= n
    v = v[1:χ, :]
    return u, s, v, n
end

function to_mps2(Θ, mps, Λ_old; χ_max=100, trim=1e-10)
    ds = [get_physical_dim(site) for site in mps]
    n = length(mps)
    half = div(n, 2)
    # the right-most Λ matrix is inherited from old mps
    Λ = get_Λ(mps[end])
    # SVD leaves the U and V matrices left-/right-orthogonalized respectively
    # right half should be right-orthogonalized
    right_sites = []
    for i=n:(-1):(half+2)
        χL, d, χR = size(Θ)
        Θ = reshape(Θ, (χL*prod(ds[1:(i-1)]), ds[i]*χR))
        u, s, v, n = decompose(Θ, χ_max, trim)
        χ = length(s)
        Γ = reshape(v, (χ, ds[i], χR))
        pushfirst!(right_sites, Site(Γ, Λ))
        Λ = Diagonal(s)
        Θ = reshape(u, (χL, prod(ds[1:(i-1)]), χ))
    end
    # the left-most Λ matrix is saved for later
    Λ_center = Λ
    # left half should be left-orthogonalized
    left_sites = []
    for i=1:half
        χL, d, χR = size(Θ)
        Θ = reshape(Θ, (χL*ds[i], prod(ds[(i+1):end])*χR))
        u, s, v, n = decompose(Θ, χ_max, trim)
        χ = length(s)
        Γ = reshape(u, (χL, ds[i], χ))
        Λ = Diagonal(s)
        push!(left_sites, Site(Γ, Λ))
        Θ = reshape(v, (χ, prod(ds[(i+1):end]), χR))
    end
    Γ = Θ
    pushfirst!(right_sites, Site(Γ, Λ_center))
    set_Γ(left_sites[1], contract_left_Γ(pinv(Λ_old), left_sites[1]))
    set_Γ(right_sites[end], contract_right_Γ(right_sites[end], pinv(get_Λ(mps[end]))))
    return MPS(vcat(left_sites, right_sites), mps.infinite)
end

function to_mps3(Θ, mps; χ_max=100, trim=1e-10)
    ds = [get_physical_dim(site) for site in mps]
    n = length(mps)
    half = div(n, 2)
    # the right-most Λ matrix is inherited from old mps
    Λ = get_Λ(mps[end])
    # SVD leaves the U and V matrices left-/right-orthogonalized respectively
    # right half should be right-orthogonalized
    right_sites = []
    for i=n:(-1):(half+2)
        χL, d, χR = size(Θ)
        Θ = reshape(Θ, (χL*prod(ds[1:(i-1)]), ds[i]*χR))
        u, s, v, n = decompose(Θ, χ_max, trim)
        χ = length(s)
        Γ = reshape(v, (χ, ds[i], χR))
        pushfirst!(right_sites, Site(Γ, Λ))
        Λ = Diagonal(s)
        Θ = reshape(u, (χL, prod(ds[1:(i-1)]), χ))
    end
    # the left-most Λ matrix is saved for later
    Λ_center = Λ
    # left half should be left-orthogonalized
    left_sites = []
    for i=1:half
        χL, d, χR = size(Θ)
        Θ = reshape(Θ, (χL*ds[i], prod(ds[(i+1):end])*χR))
        u, s, v, n = decompose(Θ, χ_max, trim)
        χ = length(s)
        Γ = reshape(u, (χL, ds[i], χ))
        Λ = Diagonal(s)
        push!(left_sites, Site(Γ, Λ))
        Θ = reshape(v, (χ, prod(ds[(i+1):end]), χR))
    end
    Γ = Θ
    pushfirst!(right_sites, Site(Γ, Λ_center))
    return MPS(vcat(left_sites, right_sites), mps.infinite)
end

function to_mps(Θ, mps::MPS; kwargs...)
    Λ_old = get_Λ(mps[end])
    mps = to_mps2(Θ, mps, Λ_old; kwargs...)
    @assert norm(mps) ≈ 1
    return mps
end

function to_mps(Θ, mps::MPS, r::UnitRange; kwargs...)
    a, b = r.start, r.stop
    Λ_old = get_left_Λ(mps, a)
    return to_mps2(Θ, mps[a:b], Λ_old; kwargs...)
end

function left_orthogonalize!(mps::MPS; χ_max=100, trim=1e-10)
    for (i, site) in enumerate(mps)
        S = get_A(site)
        χL, d, χR = size(S)
        S = reshape(S, (χL*d, χR))
        u, s, v, n = decompose(S, χ_max, trim)
        χ = length(s)
        set_Γ(mps[i], reshape(u, (χL, d, χ)))
        set_Λ(mps[i], Diagonal(s))
        if i != length(mps)
            set_Γ(mps[i+1], contract_left_Γ(v, mps[i+1]))
        else
            println("Left-orthogonalized. Final v=$v")
        end
    end
end

function right_orthogonalize!(mps::MPS; χ_max=100, trim=1e-10)
    for i=length(mps):(-1):1
        if i < 1
            S = get_Λ(mps[i-1]) * get_Γ(mps[i])
        else
            S = get_Γ(mps[i])
        end
        χL, d, χR = size(S)
        S = reshape(S, (χL, d*χR))
        u, s, v, n = decompose(S, χ_max, trim)
        χ = length(s)
        set_Γ(mps[i], reshape(v, (χ, d, χR)))
        if i != 1
            set_Λ(mps[i-1], Diagonal(s))
            set_Γ(mps[i-1], contract_right_Γ(mps[i-1], u))
        else
            println("Right-orthogonalized. Final u=$u")
        end
    end
end

function coarsegrain(mps::MPS)
    Θ = get_Γ(mps[1])
    for (i, site) in enumerate(mps[1:(end-1)])
        S = contract_left_Γ(get_Λ(mps[i]), mps[i+1])
        @tensor begin
            Θ[α, i, j, β] := Θ[α, i, γ] * S[γ, j, β]
        end
        χL, di, dj, χR = size(Θ)
        Θ = reshape(Θ, (χL, di*dj, χR))
    end
    Γ = Θ
    Λ = get_Λ(mps[end])
    return Site(Γ, Λ)
end
