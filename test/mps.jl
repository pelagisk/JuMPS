eigs_kwargs = Dict(:maxiter=>300, :tol=>0.0)
svd_kwargs = Dict(:χ_max=>50, :trim=>1e-10)

sites = []
for i=1:2
    Γ = reshape([1, 1]/sqrt(2), (1, 2, 1))
    Λ = Diagonal([1.0])
    S = Site(Γ, Λ)
    push!(sites, S)
end
mps = MPS(sites; infinite=true)

@test distance(mps, to_mps(get_Θ(mps), mps; svd_kwargs...)) < 1e-15
@test norm(mps) ≈ 1
@test norm(to_mps(get_Θ(mps), mps; svd_kwargs...)) ≈ 1

sites = []
for i=1:10
    Γ = reshape([1, 1]/sqrt(2), (1, 2, 1))
    Λ = Diagonal([1.0])
    S = Site(Γ, Λ)
    push!(sites, S)
end
mps = MPS(sites; infinite=false)

# right_orthogonalize!(mps)
# for site in mps
#     @show is_right_orthogonal(site)
# end
# left_orthogonalize!(mps)
# for i=2:length(mps)
#     site = Site(get_Γ(mps[i]), get_Λ(mps[i-1]))
#     @show is_left_orthogonal(site)
# end

for i=1:10
    Θ = rand(10, 10) + 1im * rand(10, 10)
    Θ = Θ + Θ'
    Θ /= norm(Θ)
    u, s, v, n = decompose(Θ, 1000, 1e-16)
    @test abs(norm(vec(u*Diagonal(s)*v-Θ))) < 1e-13
end
