module MPIRecoKernelAbstractionsExt 

using MPIReco, MPIReco.Adapt, MPIReco.LinearAlgebra, MPIReco.RegularizedLeastSquares, MPIReco.LinearOperatorCollection
using MPIReco.AbstractImageReconstruction
using KernelAbstractions, GPUArrays
using KernelAbstractions.Extras: @unroll
using Atomix

include("MultiPatch.jl")
include("Weighting.jl")

end