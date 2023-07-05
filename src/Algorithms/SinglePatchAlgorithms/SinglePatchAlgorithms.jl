export SinglePatchReconstructionAlgorithm, SinglePatchReconstructionParameter, SinglePatchParameters

abstract type AbstractSinglePatchReconstructionAlgorithm <: AbstractMPIReconstructionAlgorithm end
abstract type AbstractSinglePatchReconstructionParameters <: AbstractReconstructionParameters end
abstract type AbstractSinglePatchAlgorithmParameters <: AbstractMPIRecoParameters end

Base.@kwdef mutable struct SinglePatchParameters{PR<:AbstractPreProcessingParameters,
     R<:AbstractSinglePatchReconstructionParameters, PT<:AbstractPostProcessingParameters} <: AbstractSinglePatchAlgorithmParameters
  pre::PR
  reco::R
  post::PT = NoPostProcessing() 
end
  
function RecoUtils.process(algo::T, data::MPIFile, params::SinglePatchParameters) where {T<:AbstractSinglePatchReconstructionAlgorithm}
  result = process(algo, data, params.pre)
  result = process(algo, result, params.reco)
  result = process(algo, result, params.post)
  return result
end


include("SinglePatchAlgorithm.jl")
include("SinglePatchTwoStepReconstructionAlgorithm.jl")
include("SinglePatchBGEstimationAlgorithm.jl")
#include("SinglePatchTemporalRegularizationAlgorithm.jl")