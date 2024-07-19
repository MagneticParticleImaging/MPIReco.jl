module MPIRecoKernelAbstractionsExt 

using MPIReco, MPIReco.Adapt, MPIReco.LinearAlgebra
using KernelAbstractions, GPUArrays
using KernelAbstractions.Extras: @unroll

include("MultiPatch.jl")

end