# function superblock_hamiltonian(L, mpos, R)
#     χL = size(L, 2)
#     χR = size(R, 2)
#     ds = map(x -> size(x, 3), mpos)
#     function h(x::AbstractVector)::AbstractVector
#         d = ds[1]
#         x = reshape(x, (χL, d, div(length(x), d*χL)))
#         # first contract with left environment
#         @tensor begin
#             z[A, i, a, B] := L[a, C, A] * x[C, i, B]
#         end
#         # then contract with the mpos
#         for (i, mpo) in enumerate(mpos)
#             @tensor begin
#                 z[A, i, a, B] := z[A, j, b, B] * mpo[b, a, j, i]
#             end
#             # if not final iteration
#             if i != length(ds)
#                 s1, s2, s3, s4 = size(z)
#                 new_shape = (s1*s2, s3, ds[i+1], div(s4, ds[i+1]))
#                 z = reshape(z, new_shape)
#                 z = permutedims(z, (1, 3, 2, 4))
#             end
#         end
#         # contract with the right environment
#         s1, s2, s3, s4  = size(z)
#         z = reshape(z, (s1*s2, s3, s4))
#         @tensor begin
#             z[A, B] := z[A, a, C] * R[a, C, B]
#         end
#         # unfold
#         z = vec(z)
#         return z
#     end
#     dim = χL * prod(ds) * χR
#     H = LinearMap{ComplexF64}(h, dim; ismutating=false, ishermitian=true, isposdef=true)
#     return H
# end
