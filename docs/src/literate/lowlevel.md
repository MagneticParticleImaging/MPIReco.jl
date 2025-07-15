Finally, we have arrived at the low level reconstruction routine that has the signature
```julia
function reconstruction(S, u::Array; sparseTrafo = nothing,
                        lambd=0, progress=nothing, solver = "Kaczmarz",
                        weights=nothing, kargs...)
```
One can see that it requires the system matrix `S` and the measurements `u` to be
already loaded.

We note that `S` is typeless for a reason here. For a regular reconstruction one
will basically feed in an `Array{ComplexF32,2}` in here, although more precisely
it will be a `Transposed` version of that type if the `Kaczmarz` algorithm is being
used for efficiency reasons.

However, in case that matrix compression is applied `S` will be of type `SparseMatrixCSC`.
And for [Multi-Patch Reconstruction](@ref) `S` will be of type `MultiPatchOperator`. Hence,
the solvers are implemented in a very generic way and require only certain functions
to be implemented. The low level reconstruction method calls one of the solvers
from [RegularizedLeastSquares.jl](https://github.com/tknopp/RegularizedLeastSquares.jl).
