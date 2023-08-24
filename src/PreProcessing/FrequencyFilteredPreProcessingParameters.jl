export FrequencyFilteredPreProcessingParameters
Base.@kwdef struct FrequencyFilteredPreProcessingParameters{B, T<:AbstractMPIPreProcessingParameters{B}} <: AbstractMPIPreProcessingParameters{B}
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

function process(t::Type{<:AbstractMPIRecoAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParameters{NoBackgroundCorrectionParameters, <:CommonPreProcessingParameters})
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))), ignore = [:neglectBGFrames, :bgParams])
  result = getMeasurementsFD(f, bgCorrection = false; kwargs...)
  return result
end
function process(t::Type{<:AbstractMPIRecoAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParameters{InternalBackgroundCorrectionParameters, <:CommonPreProcessingParameters})
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))), ignore = [:neglectBGFrames, :bgParams])
  result = getMeasurementsFD(f; bgCorrection = true, interpolateBG = params.bgParams.interpolateBG, kwargs...)
  return result
end
function process(t::Type{<:AbstractMPIRecoAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParameters{<:ExternalBackgroundCorrection, <:CommonPreProcessingParameters})
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))), ignore = [:neglectBGFrames, :bgParams])
  bgParams = fromKwargs(FrequencyFilteredBackgroundCorrectionParameters; kwargs..., bgParams = params.bgParams)
  result = getMeasurementsFD(f, bgCorrection = false; kwargs...)
  return process(t, result, bgParams)
end