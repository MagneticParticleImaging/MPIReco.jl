export AbstractMultiPatchReconstructionAlgorithm, AbstractMultiPatchReconstructionParameters, AbstractMultiPatchAlgorithmParameters

abstract type AbstractMultiPatchReconstructionAlgorithm <: AbstractMPIReconstructionAlgorithm end
abstract type AbstractMultiPatchReconstructionParameters <: AbstractReconstructionParameters end
abstract type AbstractMultiPatchAlgorithmParameters <: AbstractMPIRecoParameters end

Base.@kwdef mutable struct MultiPatchParameters{PR<:AbstractPreProcessingParameters,
     R<:AbstractMultiPatchReconstructionParameters, PT<:AbstractPostProcessingParameters} <: AbstractMultiPatchAlgorithmParameters
  pre::PR
  reco::R
  post::PT = NoPostProcessing() 
end
  
function RecoUtils.process(algo::T, data::MPIFile, params::MultiPatchParameters) where {T<:AbstractMultiPatchReconstructionAlgorithm}
  result = process(algo, data, params.pre)
  result = process(algo, result, params.reco)
  result = process(algo, result, params.post)
  return result
end


include("MultiPatchAlgorithm.jl")
