export AbstractMultiPatchReconstructionAlgorithm, AbstractMultiPatchReconstructionParameters, AbstractMultiPatchAlgorithmParameters, MultiPatchParameters

abstract type AbstractMultiPatchReconstructionAlgorithm <: AbstractMPIRecoAlgorithm end
abstract type AbstractMultiPatchReconstructionParameters <: AbstractMPIReconstructionParameters end
abstract type AbstractMultiPatchAlgorithmParameters <: AbstractMPIRecoParameters end

export AbstractFocusFieldPositions, DefaultFocusFieldPositions, CustomFocusFieldPositions
abstract type AbstractFocusFieldPositions <: AbstractMPIRecoParameters end
struct DefaultFocusFieldPositions <: AbstractFocusFieldPositions end
positions(ffPos::DefaultFocusFieldPositions) = nothing
@parameter struct CustomFocusFieldPositions{T<:AbstractArray} <: AbstractFocusFieldPositions
  positions::T
end
positions(ffPos::CustomFocusFieldPositions) = ffPos.positions

@chain mutable struct MultiPatchParameters{PR<:AbstractMPIPreProcessingParameters,
     R<:AbstractMultiPatchReconstructionParameters, PT<:AbstractMPIPostProcessingParameters} <: AbstractMultiPatchAlgorithmParameters
  pre::Union{PR, AbstractUtilityReconstructionParameters{PR}}
  reco::R
  post::PT = NoPostProcessing() 
end
  
function (params::MultiPatchParameters)(algo::T, data::AbstractArray, args...) where {T<:AbstractMultiPatchReconstructionAlgorithm}
  throw(ArgumentError("MultiPatchAlgorithms are not defined for the given arguments, expected <: MPIFile, found $(typeof(data))"))
end



include("MultiPatchAlgorithm.jl")
include("MultiPatchPeriodicMotion.jl")