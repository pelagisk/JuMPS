"
Here goes code that defines a model, which should
- be practical from a user perspective
- generate MPO's
- generate local Hamiltonians

Will probably use Crosswhite's finite state automata.
"

struct Model
    mps
    terms
end

function find_terms(model, i)
    terms = Dict()
    # TODO use filter
    for (key, value) in pairs(model.terms)
        if isinfinite(model.mps)
            # @show key
            key = mod.(key .- 1, length(model.mps)) .+ 1
            # @show key
        end
        if (i in key)
            terms[key] = value
        end
    end
    return terms
end

function mpos(model::Model)
    Ms = []
    for (i, site) in enumerate(model.mps)
        # TODO check this here
        @show i
        # find all the relevant terms
        # TODO get position in terms here as well?
        terms = find_terms(model, i)
        # longest tuple of site indices sets the MPO virtual dimension?
        if length(terms) > 0
            D = maximum(length.(keys(terms))) + 1
            d = get_physical_dim(site)
            mpo = zeros(ComplexF64, D, D, d, d)
            for (positions, v) in pairs(terms)
                cnumber = v[1]
                opkeys = v[2:end]
                @show opkeys
                j = findfirst(x -> x == i, positions)
                x = D - length(opkeys) + j
                y = j
                println("placing $(opkeys[j]) at position $x, $y")
                mpo[x, y, :, :] = get_operator(site, opkeys[j])
                mpo
            end
            push!(Ms, mpo)
        end
    end
    return Ms
end
