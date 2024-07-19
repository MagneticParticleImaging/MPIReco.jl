module MPIRecoKernelAbstractionsExt 

using MPIReco, MPIReco.Adapt, MPIReco.LinearAlgebra
using KernelAbstractions, GPUArrays
using KernelAbstractions.Extras: @unroll
using KernelAbstractions.Atomix

include("MultiPatch.jl")

end