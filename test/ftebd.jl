# at criticality
J = 1.0
h = J
n_sites = 10

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
    Γ = reshape([1, -1]/sqrt(2), (1, 2, 1))
    Λ = Diagonal([1.0])
    S = Site(Γ, Λ)
    push!(sites, S)
end
mps = MPS(sites; infinite=false)

svd_kwargs = Dict(:χ_max=>50, :trim=>1e-10)

δt = 1e-1
U_local = exp(-δt * H_local)

etol = 1e-8
e, mps = ftebd(mps, U_local, H_local; verbose=2, etol=etol, svd_kwargs=svd_kwargs)
@test abs(e - ising_exact_energy_critical_finite(J, n_sites)) < 1e-3

# TODO extremely ineffective - why?
