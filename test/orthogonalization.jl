J = 0.1
h = 1.0
n_sites = 20

function ising_exact_energy_critical_finite(J, n_sites)
    return J*(1 - 1/sin(π/(2*(2*n_sites+1))))
end

σx = [0 1; 1 0]
σz = [1 0; 0 -1]
o = [0 0; 0 0]
i = [1 0; 0 1]

mpo = zeros(3, 3, 2, 2)
mpo[1, 1, :, :] = i
mpo[3, 3, :, :] = i
mpo[3, 1, :, :] = h * σx
mpo[2, 1, :, :] = -J * σz
mpo[3, 2, :, :] = σz
mpos = [mpo for i=1:n_sites]

H_local = -J * kron(σz, σz) + h/2 * (kron(σx, i) + kron(i, σx))

sites = []
for i=1:n_sites
    Γ = reshape([1, 1]/sqrt(2), (1, 2, 1))
    Λ = Diagonal([1.0])
    S = Site(Γ, Λ)
    push!(sites, S)
end
mps = MPS(sites; infinite=false)

eigs_kwargs = Dict(:maxiter=>300, :tol=>0.0)
svd_kwargs = Dict(:χ_max=>50, :trim=>1e-10)

etol = 1e-5
e, mps = fdmrg(mps, mpos, H_local; verbose=1, etol=etol, eigs_kwargs=eigs_kwargs, svd_kwargs=svd_kwargs)

left_orthogonalize!(mps)
for site in mps
    @test is_left_orthogonal(site) == true
end

right_orthogonalize!(mps)
for i=2:length(mps)
    site = Site(get_Γ(mps[i]), get_Λ(mps[i-1]))
    @test is_right_orthogonal(site) == true
end
