export CommonPreProcessingParameters
# TODO Averaging:
# numAveraging always
# numPeriodAverages | averagePeriodsPerPatch, potentially as small structs?
Base.@kwdef struct CommonPreProcessingParameters{T<:AbstractBackgroundCorrectionParameters} <: AbstractPreProcessingParameters{T}
  bgParams::T = NoBackgroundCorrectionParameters()
  neglectBGFrames::Bool = true
  numAverages::Int64 = 1
  numPeriodAverages::Int64 = 1
  numPeriodGrouping::Int64 = 1
  frames::Union{Nothing, UnitRange{Int64}, Vector{Int64}} = nothing
  spectralLeakageCorrection::Bool = false
  loadasreal::Bool = false
end

function process(t::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::CommonPreProcessingParameters{NoBackgroundCorrectionParameters})
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))), ignore = [:neglectBGFrames, :bgParams])
  result = getMeasurementsFD(f; bgCorrection = false, kwargs...)
  return result
end
function process(t::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::CommonPreProcessingParameters{InternalBackgroundCorrectionParameters})
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))), ignore = [:neglectBGFrames, :bgParams])
  result = getMeasurementsFD(f; bgCorrection = true, interpolateBG = params.bgCorrection.interpolateBG, kwargs...)
  return result
end
function process(t::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::CommonPreProcessingParameters{<:ExternalBackgroundCorrection})
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))), ignore = [:neglectBGFrames, :bgParams])
  result = getMeasurementsFD(f, bgCorrection = false, kwargs...)
  bgParams = fromKwargs(ExternalPreProcessedBackgroundCorrectionParameters; kwargs..., bgParams = params.bgParams)
  return process(t, result, bgParams)
end