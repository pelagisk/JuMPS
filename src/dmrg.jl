function groundstate(H, χL, d, χR; kwargs...)
    e, v = eigs(H; nev=1, which=:SR, kwargs...)
    e = real(e[1])
    v /= norm(v)
    Θ = reshape(v, (χL, d, χR))
    return e, Θ
end

function groundstate(H, Θ; kwargs...)
    e, v = eigs(H; nev=1, which=:SR, v0=vec(Θ), kwargs...)
    e = real(e[1])
    v /= norm(v)
    Θ = reshape(v, size(Θ))
    return e, Θ
end

function superblock_hamiltonian(L, mpos, R)
    χL = size(L, 2)
    χR = size(R, 2)
    # contract the MPOs
    mpo = mpos[1]
    for i=2:length(mpos)
        @tensor begin
            mpo[A, B, i, j, k, l] := mpo[A, C, i, k] * mpos[i][C, B, j, l]
        end
        DA, DB, di, dj, dk, dl = size(mpo)
        mpo = reshape(mpo, (DA, DB, di*dj, dk*dl))
    end
    d = size(mpo)[end]
    # define the linear operator
    function h(x::AbstractVector)::AbstractVector
        Θ = reshape(x, (χL, d, χR))
        @tensor begin
            Θ[α, a, β] := L[A, γ, α] * Θ[γ, b, δ] * R[B, δ, β] * mpo[A, B, b, a]
        end
        return vec(Θ)
    end
    dim = χL * d * χR
    H = LinearMap{ComplexF64}(h, dim; ismutating=false, ishermitian=true, isposdef=true)
    return H
end

function update_left_environment!(L, mpos, mps, r::UnitRange)
    a, b = r.start, r.stop
    n = b-a+1
    half = div(n, 2)
    for i=a:(a+half-1)
        S = contract_left_Γ(get_left_Λ(mps, i), mps[i])
        @tensor begin
            L[A, α, β] := L[B, γ, δ] * S[γ, a_, α] * mpos[i][B, A, a_, b] * conj(S)[δ, b, β]
        end
    end
    return L
end

function update_left_environment!(L, mpos, mps)
    return update_left_environment!(L, mpos, mps, 1:length(mps))
end

function update_right_environment!(R, mpos, mps, r::UnitRange)
    a, b = r.start, r.stop
    n = b-a+1
    half = div(n, 2)
    for i=b:(-1):(b-half+1)
        S = get_A(mps[i])
        @tensor begin
            R[A, α, β] := S[α, a, γ] * mpos[i][A, B, a, b] * conj(S)[β, b, δ] * R[B, γ, δ]
        end
    end
    return R
end

function update_right_environment!(R, mpos, mps)
    return update_right_environment!(R, mpos, mps, 1:length(mps))
end

"
TODO
- iteration measurements
- return DMRG object or dict?
- start from general input MPS?
- connect to model definition
"

function idmrg(mps, mpos, H_local; verbose=false, maxiter=1000, etol=1e-12, stol=1e-8, eigs_kwargs=Dict(maxiter=300, tol=0.0), svd_kwargs=Dict(χ_max=100, trim=1e-10))

    L = zeros(ComplexF64, size(mpos[1], 1), 1, 1)
    L[end, 1, 1] = 1
    R = zeros(ComplexF64, size(mpos[end], 1), 1, 1)
    R[1, 1, 1] = 1

    entanglements = [0.0]; δs = 1.0
    energies = [0.0]; δe = 1.0
    i = 1
    log_print(verbose,
        "i\tχ\ts\t\t\te",
        "i\tχ\ts\t\t\tδs\t\t\te\t\t\tδe",
    )
    while (i < maxiter) && ((δs > stol) || (δe > etol))
        # @assert distance(mps, to_mps(get_Θ(mps), mps; svd_kwargs...)) < 1e-3
        Θ = get_Θ(mps)
        sh = superblock_hamiltonian(L, mpos, R)
        e, Θ = groundstate(sh, Θ; eigs_kwargs...)
        mps = to_mps(Θ, mps; svd_kwargs...)
        e = expectation(mps, H_local)
        push!(energies, e)
        δe = abs(e-energies[i])
        χ = maximum(get_bond_dims(mps))
        s = maximum(get_entanglements(mps))
        push!(entanglements, s)
        δs = abs(s-entanglements[i])
        log_print(verbose,
            (@sprintf "%d\t%d\t%0.16e\t%0.16e" i χ s e),
            (@sprintf "%d\t%d\t%0.16e\t%0.16e\t%0.16e\t%0.16e" i χ s δs e δe),
        )
        L = update_left_environment!(L, mpos, mps)
        R = update_right_environment!(R, mpos, mps)
        mps = shuffle!(mps)
        i += 1
    end

    return energies[end], mps
end

function fdmrg(mps, mpos, H_local; verbose=false, focus=2, maxiter=1000, etol=1e-12, stol=1e-8, eigs_kwargs=Dict(maxiter=300, tol=0.0), svd_kwargs=Dict(χ_max=100, trim=1e-10))

    leftsweep = collect(partition(1:length(mps), focus, 1))
    rightsweep = reverse(leftsweep[1:(end-1)])
    sweep = vcat(leftsweep, rightsweep)

    Ls = OrderedDict()
    Rs = OrderedDict()
    for (a, b) in leftsweep
        L = zeros(ComplexF64, size(mpos[1], 1), 1, 1)
        L[end, 1, 1] = 1
        Ls[a] = L
        R = zeros(ComplexF64, size(mpos[end], 1), 1, 1)
        R[1, 1, 1] = 1
        Rs[b] = R
    end
    entanglements = [0.0]; δs = 1.0
    energies = [0.0]; δe = 1.0
    i = 1
    while (i < maxiter) && ((δs > stol) || (δe > etol))
        log_print(verbose, "--- sweep $i ---")
        log_print(verbose,
            "(a, b)\tχ\ts\t\t\te",
            "(a, b)\tχ\ts\t\t\tδs\t\t\te\t\t\tδe",
        )
        for (j, (a, b)) in enumerate(sweep)
            Θ = get_Θ(mps, a:b)
            sh = superblock_hamiltonian(Ls[a], mpos[a:b], Rs[b])
            e, Θ = groundstate(sh, Θ; eigs_kwargs...)
            mps[a:b] = to_mps(Θ, mps, a:b; svd_kwargs...)
            # e = expectation(mps, H_local)
            push!(energies, e)
            δe = e-energies[i]
            χ = maximum(get_bond_dims(mps))
            s = maximum(get_entanglements(mps))
            push!(entanglements, s)
            δs = s-entanglements[i]
            log_print(verbose,
                (@sprintf "(%d, %d)\t%d\t%0.16e\t%0.16e" a b χ s e),
                (@sprintf "(%d, %d)\t%d\t%0.16e\t%0.16e\t%0.16e\t%0.16e" a b χ s δs e δe),
            )
            # if we are sweeping right
            if (j < length(leftsweep))
                Ls[a+1] = update_left_environment!(Ls[a], mpos, mps, a:b)
            else  # or left
                Rs[b-1] = update_right_environment!(Rs[b], mpos, mps, a:b)
            end
        end
        i += 1
    end
    return energies[end], mps
end
