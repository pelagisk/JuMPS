J = 0.8
h = 1.0

function ising_exact_energy(J, h)
    J = -abs(J)
    E_exact, error = quadgk(k -> J*sqrt(1 + (h/J)^2 - 2*(h/J)*cos(k))/pi, 0, pi)
    return E_exact
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
mpos = [mpo, mpo]

H_local = -J * kron(σz, σz) + h/2 * (kron(σx, i) + kron(i, σx))

sites = []
for i=1:2
    Γ = reshape([1, 1]/sqrt(2), (1, 2, 1))
    Λ = Diagonal([1.0])
    S = Site(Γ, Λ)
    push!(sites, S)
end
mps = MPS(sites; infinite=true)

eigs_kwargs = Dict(:maxiter=>300, :tol=>0.0)
svd_kwargs = Dict(:χ_max=>50, :trim=>1e-10)

etol = 1e-9
e, mps = idmrg(mps, mpos, H_local; verbose=1, etol=etol, eigs_kwargs=eigs_kwargs, svd_kwargs=svd_kwargs)
@test abs(e - ising_exact_energy(J, h)) < etol
@test abs(norm_and_correlation_length(mps)[1] - 1.0) < 1e-12

# TODO test canonicalize?

canonicalize!(mps)
