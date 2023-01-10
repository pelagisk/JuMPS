module Tst

using Test, DataStructures, LinearAlgebra, QuadGK

include("../src/JuMPS.jl")
using .JuMPS

include("mps.jl")
include("transfer.jl")
# include("itebd.jl")
# include("ftebd.jl")
# include("idmrg.jl")
# include("fdmrg.jl")
# include("orthogonalization.jl")
include("model.jl")

# TODO tests for using unit cell of other length!

end
