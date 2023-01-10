# AKLT state
sites = []
for i=1:2
    Γ = zeros(2, 3, 2)
    Γ[:, 1, :] = sqrt(4/3) * [0 1; 0 0]
    Γ[:, 2, :] = -sqrt(2/3) * [1 0; 0 -1]
    Γ[:, 3, :] = -sqrt(4/3) * [0 0; 1 0]
    Λ = Diagonal([1, 1]/sqrt(2))
    S = Site(Γ, Λ)
    push!(sites, S)
end
mps = MPS(sites; infinite=true)

# Product state
sites = []
for i=1:2
    Γ = reshape([1, 1]/sqrt(2), (1, 2, 1))
    Λ = Diagonal([1.0])
    S = Site(Γ, Λ)
    push!(sites, S)
end
mps = MPS(sites; infinite=true)

canonicalize!(mps)
