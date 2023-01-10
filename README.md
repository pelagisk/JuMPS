# JuMPS

Yet another matrix product state (MPS) implementation, written in Julia.

# Usage

Install `julia`, version `> 1.0` and the packages `Arpack, TensorOperations, LinearMaps`. Then, in the Julia REPL:

```julia
include("../path/to/JuMPS.jl")
using .JuMPS
```

The first run will be slow, since Julia compiles the functions.

# Roadmap

- [ ] model -> MPOs
- [ ] exact diagonalization
- [x] MPS
  - [x] left-/right-orthogonalization
  - [x] use svds when possible
  - [ ] factorize: multiple-dispatch used for different MPS types
  - [ ] symmetric MPS
- [x] infinite DMRG
- [x] finite DMRG
- [x] infinite TEBD
- [x] finite TEBD
- [ ] infinite TDVP
- [ ] finite TDVP
- [x] transfer matrices
  - [x] canonicalization
  - [x] correlation length
  - [ ] moments


# References

Under construction!
