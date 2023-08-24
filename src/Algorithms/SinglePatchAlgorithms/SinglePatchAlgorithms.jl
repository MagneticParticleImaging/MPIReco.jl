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
  
function process(algo::T, data::MPIFile, params::SinglePatchParameters) where {T<:AbstractSinglePatchReconstructionAlgorithm}
  result = process(algo, data, params.pre)
  result = process(algo, result, params.reco)
  result = process(algo, result, params.post)
  return result
end


include("SinglePatchAlgorithm.jl")
include("SinglePatchTwoStepReconstructionAlgorithm.jl")
include("SinglePatchBGEstimationAlgorithm.jl")
#include("SinglePatchTemporalRegularizationAlgorithm.jl")