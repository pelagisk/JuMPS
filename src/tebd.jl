function propagate(Θ, U_local)
    @tensor begin
        Θ[α, a, β] :=  U_local[a, b] * Θ[α, b, β]
    end
    return Θ
end

function itebd(mps, U_local, H_local; verbose=false, maxiter=1E8, etol=1e-12, stol=1e-8, svd_kwargs=Dict(χ_max=100, trim=1e-10))
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
        Θ = propagate(Θ, U_local)
        Θ /= norm(vec(Θ))
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
        mps = shuffle!(mps)
        i += 1
    end

    return energies[end], mps
end

function ftebd(mps, U_local, H_local; verbose=false, focus=2, maxiter=1E8, etol=1e-12, stol=1e-8, svd_kwargs=Dict(χ_max=100, trim=1e-10))
    leftsweep = collect(partition(1:length(mps), focus, 1))
    rightsweep = reverse(leftsweep[1:(end-1)])
    sweep = vcat(leftsweep, rightsweep)
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
            Θ = propagate(Θ, U_local)
            Θ /= norm(vec(Θ))
            mps[a:b] = to_mps(Θ, mps, a:b; svd_kwargs...)
            # e = expectation(mps, H_local)
            e = energies[end]
            # push!(energies, e)
            # δe = e-energies[i]
            χ = maximum(get_bond_dims(mps))
            s = maximum(get_entanglements(mps))
            push!(entanglements, s)
            δs = s-entanglements[i]
            log_print(verbose,
                (@sprintf "(%d, %d)\t%d\t%0.16e\t%0.16e" a b χ s e),
                (@sprintf "(%d, %d)\t%d\t%0.16e\t%0.16e\t%0.16e\t%0.16e" a b χ s δs e δe),
            )
        end
        i += 1
    end
    return energies[end], mps
end
