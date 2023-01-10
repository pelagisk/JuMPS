"
Here goes code that defines transfer matrices and operations with them.

This code will be used for:
- correlation lengths
- moments of observables, such as variance of the Hamiltonian
- normalization
- calculation of observables
"

mutable struct Transfer
    Θ1::Array{ComplexF64, 3}
    operator::Array{ComplexF64, 2}
    Θ2::Array{ComplexF64, 3}
    identity_operator::Bool
    conjugate::Bool
end

function Transfer(mps::MPS; right=true)
    if right == true
        Θ = get_A(mps[1])
        for (i, site) in enumerate(mps[2:end])
            @tensor begin
                Θ[α, i, j, β] := Θ[α, i, γ] * get_A(site)[γ, j, β]
            end
            χL, di, dj, χR = size(Θ)
            Θ = reshape(Θ, (χL, di*dj, χR))
        end
    else  # left transfer
        Θ = if isinfinite(mps)
            get_Λ(mps[end])
        else
            Diagonal([1.0])
        end
        Θ = contract_left_Γ(Θ, mps[1])
        for (i, site) in enumerate(mps[1:(end-1)])
            S = contract_left_Γ(get_Λ(mps[i]), get_Γ(mps[i+1]))
            @tensor begin
                Θ[α, i, j, β] := Θ[α, i, γ] * S[γ, j, β]
            end
            χL, di, dj, χR = size(Θ)
            Θ = reshape(Θ, (χL, di*dj, χR))
        end
    end
    return Transfer(Θ, zeros(1, 1), zeros(1, 1, 1), true, true)
end

function Transfer(mps::MPS, r::UnitRange; right=true)
    a, b = r.start, r.stop
    if right == true
        Θ = get_A(mps[a])
        for (i, site) in enumerate(mps[(a+1):b])
            @tensor begin
                Θ[α, i, j, β] := Θ[α, i, γ] * get_A(site)[γ, j, β]
            end
            χL, di, dj, χR = size(Θ)
            Θ = reshape(Θ, (χL, di*dj, χR))
        end
    else  # left transfer
        Θ = get_left_Λ(mps, a)
        Θ = contract_left_Γ(Θ, mps[1])
        for (i, site) in enumerate(mps[1:(end-1)])
            S = contract_left_Γ(get_Λ(mps[i]), get_Γ(mps[i+1]))
            @tensor begin
                Θ[α, i, j, β] := Θ[α, i, γ] * S[γ, j, β]
            end
            χL, di, dj, χR = size(Θ)
            Θ = reshape(Θ, (χL, di*dj, χR))
        end
    end
    return Transfer(Θ, zeros(1, 1), zeros(1, 1, 1), true, true)
end

function Transfer(site::Site; left=false)
    if left == false
        Θ = get_A(site)
    else  # left transfer
        Θ = contract_left_Γ(get_Λ(site), site)
    end
    return Transfer(Θ, zeros(1, 1), zeros(1, 1, 1), true, true)
end

# Transfer(mps::MPS) = Transfer(get_Θ(mps), zeros(1, 1), zeros(1, 1, 1), true, true)
# Transfer(mps::MPS, O) = Transfer(get_Θ(mps), O, zeros(1, 1, 1), false, true)

function Base.size(T::Transfer)
    return if T.conjugate == true
        size(T.Θ1, 1)^2, size(T.Θ1, 3)^2
    else
        size(T.Θ1, 1) * size(T.Θ2, 1), size(T.Θ1, 3) * size(T.Θ2, 3)
    end
end


function operator_representation(T::Transfer)
    S = T.Θ1
    # sandwich with Θ1 if Transfer is with its conjugate
    Z = (T.conjugate == true) ? T.Θ1 : T.Θ2
    @assert size(S) == size(Z)
    χL, d, χR = size(S)
    if χL < χR
        println("right larger than left")
        # TODO this is not the right way to extend, judging from the ON tests
        S_ = zeros(ComplexF64, χR, d, χR)
        S_[1:χL, :, :] = S
        S = S_
        Z_ = zeros(ComplexF64, χR, d, χR)
        Z_[1:χL, :, :] = Z
        Z = Z_
    elseif χL > χR
        println("left larger than right")
        S = S[1:χR, :, :]
        Z = Z[1:χR, :, :]
    end

    # @assert χL ≡ χR
    if T.identity_operator == true
        function ti(x::AbstractVector)::AbstractVector
            # @assert χR*χR == length(x)
            x = reshape(x, (χR, χR))
            @tensor begin
                y[α, β] := S[α, i, γ] * conj(Z[β, i, δ]) * x[γ, δ]
            end
            y = vec(y)
            return y
        end
        dim = χR^2
        operator = LinearMap{ComplexF64}(ti, dim; ismutating=false, ishermitian=true, isposdef=true)
        return operator
    else
        function t(x::AbstractVector)::AbstractVector
            # @assert χR*χR == length(x)
            x = reshape(x, (χR, χR))
            @tensor begin
                y[α, β] := S[α, i, γ] * conj(Z[β, j, δ]) * T.operator[i, j] * x[γ, δ]
            end
            y = vec(y)
            return y
        end
        dim = χR^2
        operator = LinearMap{ComplexF64}(t, dim; ismutating=false, ishermitian=false, isposdef=true)
        return operator
    end
end

function LinearAlgebra.tr(T::Transfer)::Float64
    S = T.Θ1
    # sandwich with Θ1 if Transfer is with its conjugate
    Z = (T.conjugate == true) ? T.Θ1 : T.Θ2
    # TODO good idea to keep it like this at all?
    χLs, d, χRs = size(S)
    χLz, d, χRz = size(Z)
    χ1 = min(χLs, χLz)
    χ2 = min(χRs, χRz)
    t = 0
    if T.identity_operator ≡ true
        for i=1:d
            t += LinearAlgebra.tr(S[1:χ1, i, 1:χ2] * Z[1:χ1, i, 1:χ2]')
        end
    else
        for i=1:d
            for j=1:d
                t += LinearAlgebra.tr(S[1:χ1, i, 1:χ2] * Z[1:χ1, j, 1:χ2]') * T.operator[i, j]
            end
        end
    end
    return real(t)
end

function eigen(T::Transfer; left=false, nev=1, which=:LM, kwargs...)
    if left == true
        T.Θ1 = permutedims(T.Θ1, (3, 2, 1))
        T.Θ2 = permutedims(T.Θ2, (3, 2, 1))
    end
    t = operator_representation(T)
    if size(t, 1) > 10
        e, v = eigs(t; nev=nev, which=which, kwargs...)
    else
        M = Matrix(t)
        e, v = eigen(M)
        # TODO sort after 'which' keyword
        e = e[1:nev]
        v = v[:, 1:nev]
    end
    if left == true
        T.Θ1 = permutedims(T.Θ1, (3, 2, 1))
        T.Θ2 = permutedims(T.Θ2, (3, 2, 1))
    end
    vs = []
    for j=1:size(v, 2)
        χ = Int(sqrt(length(v[:, j])))
        push!(vs, reshape(v[:, j], (χ, χ)))
    end
    return e, vs
end

function norm_and_correlation_length(mps::MPS)
    T = Transfer(mps)
    e, _ = eigen(T; nev=2, which=:LM)
    e = sort(real.(e))
    𝒩 = abs(e[2])
    ξ = log(1/abs(e[1]))
    return 𝒩, ξ
end

function make_hermitian(m)
    m = m + m'
    m /= norm(vec(m))
    if isposdef(m) == false
        m .*= -1
    end
    return m
end

function canonicalize!(mps::MPS)
    @assert isinfinite(mps) == true
    site = coarsegrain(mps)
    Tr = Transfer(site)
    Tl = Transfer(site; left=true)
    er, vr = eigen(Tr; nev=1, which=:LM); vr = vr[1]
    # these should be the left eigenvectors
    el, vlt = eigen(Tl; left=true, nev=1, which=:LM); vlt = vlt[1]
    # @show er, el
    X = cholesky(make_hermitian(vr)).L
    # @assert norm(vec(X*X' .- make_hermitian(vr))) < 1e-15
    Yt = cholesky(make_hermitian(vlt)).L
    # @assert norm(vec(Yt*Yt' .- make_hermitian(vlt))) < 1e-15
    W = Yt * get_Λ(site) * X
    F = if size(W, 1) > 20
        svds(W)[1]
    else
        svd(W)
    end
    U, λ, V = F.U, F.S, F.Vt
    Λ_new = Diagonal(λ)
    Γ_new = contract_left_right_Γ(V * inv(X), site, inv(Yt) * U)
    site = Site(Γ_new, Λ_new)
    mps = to_mps3(Γ_new, mps)
    set_Λ(mps[end], Λ_new)
    return mps
end

function is_left_orthogonal(site::Site; tol=1e-12)
    Tr = Transfer(site)
    er, vr = eigen(Tr; nev=1, which=:LM); vr = vr[1]
    χ = size(vr, 1)
    return norm(vec(abs.(vr) .- (1/sqrt(χ))*Matrix{ComplexF64}(I, χ, χ))) < tol
end

function is_right_orthogonal(site::Site; tol=1e-12)
    # note, here "site" refers to Λ_{i-1}, Γ_{i} pair
    Tl = Transfer(site; left=true)
    el, vlt = eigen(Tl; left=true, nev=1, which=:LM); vlt = vlt[1]
    χ = size(vlt, 1)
    return norm(vec(abs.(vlt) .- (1/sqrt(χ))*Matrix{ComplexF64}(I, χ, χ))) < tol
end
