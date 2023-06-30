# Not yet happy with field replication, but nested structs with parametric types is also quite ugly
export FrequencyFilteredPreProcessingParameters
Base.@kwdef struct FrequencyFilteredPreProcessingParameters{T<:AbstractBackgroundCorrectionParameters} <: AbstractPreProcessingParameters
  frequencies::Vector{Int64}
  neglectBGFrames::Bool = true
  bgCorrection::T = NoBackgroundCorrectionParameters()
  numAverages::Int64 = 1
  numPeriodAverages::Int64 = 1
  numPeriodGrouping::Int64 = 1
  frames::Union{Nothing, UnitRange{Int64}, Vector{Int64}} = nothing
  spectralLeakageCorrection::Bool = false
  loadasreal::Bool = false
end

function RecoUtils.process(t::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParameters{NoBackgroundCorrectionParameters})
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))), ignore = [:neglectBGFrames, :bgCorrection])
  result = getMeasurementsFD(f, bgCorrection = false; kwargs...)
  return result
end
function RecoUtils.process(t::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParameters{InternalBackgroundCorrectionParameters})
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))), ignore = [:neglectBGFrames, :bgCorrection])
  result = getMeasurementsFD(f; bgCorrection = true, interpolateBG = params.bgCorrection.interpolateBG, kwargs...)
  return result
end
function RecoUtils.process(t::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParameters{<:ExternalBackgroundCorrection})
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))), ignore = [:neglectBGFrames, :bgCorrection])
  bgParams = fromKwargs(FrequencyFilteredBackgroundCorrectionParameters; kwargs..., bgParams = params.bgCorrection)
  result = getMeasurementsFD(f, bgCorrection = false; kwargs...)
  return process(t, result, bgParams)
end