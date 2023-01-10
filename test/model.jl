operators = Dict(
    :σx => [0 1; 1 0],
    :σy => [0 -1im; 1im 0],
    :σz => [1 0; 0 -1],
    # :O = [0 0; 0 0],
    # :I = [1 0; 0 1]
)

SpinHalfSite = SiteType("Spin 1/2", 2, operators)

function boson_site(n)
    B = diagm(1 => sqrt.(1:n))
    N = diagm(0 => 0:n)
    operators = Dict(
        :B => B,
        :Bh => B',
        :N => N,
        :X => (B .+ B')./sqrt(2),
        :P => -1im .* (B .- B')./sqrt(2),
    )
    return SiteType("Boson (trunc=$n)", n+1, operators)
end

println("TFIM")

sites = []
for i=1:2
    Γ = reshape([1, 1]/sqrt(2), (1, 2, 1))
    Λ = Diagonal([1.0])
    S = Site(Γ, Λ, SpinHalfSite)
    push!(sites, S)
end
mps = MPS(sites; infinite=true)

h = 0.1
J = 0.5

terms = Dict()
for i=1:2
    terms[i] = (h, :σz)
    terms[(i, i+1)] = (-J, :σx, :σx)
end

model = Model(mps, terms)

Ms = mpos(model)
for (i, M) in enumerate(Ms)
    @show i
    @show size(M)
end


println("Unit cell 2, spin-boson model")

BosonSite = boson_site(5)

n_sites = 4

sites = []
for i=1:n_sites
    if i % 2 == 1
        Γ = reshape([1, 1]/sqrt(2), (1, 2, 1))
        Λ = Diagonal([1.0])
        S = Site(Γ, Λ, SpinHalfSite)
    else
        ψ = zeros(ComplexF64, 5+1)
        ψ[1] = 1.0
        Γ = reshape(ψ, (1, 5+1, 1))
        Λ = Diagonal([1.0])
        S = Site(Γ, Λ, BosonSite)
    end
    push!(sites, S)
end
mps = MPS(sites; infinite=true)

h = 0.1
ω = 1.0
J = 0.5

terms = Dict()
for i=1:n_sites
    if i % 2 == 1
        terms[i] = (h, :σz)
        terms[(i, i+1)] = (-J, :σx, :X)
    else
        terms[i] = (ω, :N)
        terms[(i, i+1)] = (-J, :X, :σx)
    end
end

model = Model(mps, terms)

Ms = mpos(model)
for (i, M) in enumerate(Ms)
    @show i
    @show size(M)
end
