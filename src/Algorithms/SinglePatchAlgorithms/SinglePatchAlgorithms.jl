export SinglePatchReconstructionAlgorithm, SinglePatchReconstructionParameter, SinglePatchParameters

abstract type AbstractSinglePatchReconstructionAlgorithm <: AbstractMPIRecoAlgorithm end
abstract type AbstractSinglePatchReconstructionParameters <: AbstractMPIReconstructionParameters end
abstract type AbstractSinglePatchAlgorithmParameters <: AbstractMPIRecoParameters end

Base.@kwdef mutable struct SinglePatchParameters{PR<:AbstractMPIPreProcessingParameters,
     R<:AbstractSinglePatchReconstructionParameters, PT<:AbstractMPIPostProcessingParameters} <: AbstractSinglePatchAlgorithmParameters
  pre::PR
  reco::R
  post::PT = NoPostProcessing() 
end
  
function process(algo::T, params::SinglePatchParameters, data::MPIFile, frequencies::Union{Vector{CartesianIndex{2}}, Nothing} = nothing) where {T<:AbstractSinglePatchReconstructionAlgorithm}
  result = process(algo, params.pre, data, frequencies)
  result = process(algo, params.reco, result)
  result = process(algo, params.post, result)
  return result
end


include("SinglePatchAlgorithm.jl")
include("SinglePatchTwoStepReconstructionAlgorithm.jl")
include("SinglePatchBGEstimationAlgorithm.jl")
#include("SinglePatchTemporalRegularizationAlgorithm.jl")