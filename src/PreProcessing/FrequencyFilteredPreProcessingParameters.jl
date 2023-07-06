export FrequencyFilteredPreProcessingParameters
Base.@kwdef struct FrequencyFilteredPreProcessingParameters{T<:AbstractPreProcessingParameters} <: AbstractPreProcessingParameters
  frequencies::Vector{Int64}
  pre::T
end

function Base.getproperty(value::FrequencyFilteredPreProcessingParameters, name::Symbol)
  if name == :frequencies || name == :pre
    return getfield(value, name)
  else
    return getproperty(getfield(value, :pre), name)
  end
end

function RecoUtils.process(t::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParameters{CommonPreProcessingParameters{NoBackgroundCorrectionParameters}})
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))), ignore = [:neglectBGFrames, :bgCorrection])
  result = getMeasurementsFD(f, bgCorrection = false; kwargs...)
  return result
end
function RecoUtils.process(t::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParameters{CommonPreProcessingParameters{InternalBackgroundCorrectionParameters}})
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))), ignore = [:neglectBGFrames, :bgCorrection])
  result = getMeasurementsFD(f; bgCorrection = true, interpolateBG = params.bgCorrection.interpolateBG, kwargs...)
  return result
end
function RecoUtils.process(t::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParameters{CommonPreProcessingParameters{<:ExternalBackgroundCorrection}})
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))), ignore = [:neglectBGFrames, :bgCorrection])
  bgParams = fromKwargs(FrequencyFilteredBackgroundCorrectionParameters; kwargs..., bgParams = params.bgCorrection)
  result = getMeasurementsFD(f, bgCorrection = false; kwargs...)
  return process(t, result, bgParams)
end