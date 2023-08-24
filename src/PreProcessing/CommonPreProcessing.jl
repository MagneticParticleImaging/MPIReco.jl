export CommonPreProcessingParameters
# TODO Averaging:
# numAveraging always
# numPeriodAverages | averagePeriodsPerPatch, potentially as small structs?
"""
Parameters for retrieving MPI measurement data in frequency domain
"""
Base.@kwdef struct CommonPreProcessingParameters{T<:AbstractMPIBackgroundCorrectionParameters} <: AbstractMPIPreProcessingParameters{T}
  "Type of background correction"
  bgParams::T = NoBackgroundCorrectionParameters()
  "Frames to be retrieved from the data"
  frames::Union{Nothing, UnitRange{Int64}, Vector{Int64}} = nothing
  "If background frames should be neglected in frame index selection"
  neglectBGFrames::Bool = true
  "Number of frames that should be averaged together"
  numAverages::Int64 = 1
  "Number of periods within a frame that should be averaged together"
  numPeriodAverages::Int64 = 1
  "Number of periods within a frame that should be grouped together, i.e. interpreted as one period"
  numPeriodGrouping::Int64 = 1
  "If spectral leakage should be corrected"
  spectralLeakageCorrection::Bool = false
  "If measurement data should be loaded as two real numbers instead of a complex one"
  loadasreal::Bool = false
end

function process(t::Type{<:AbstractMPIRecoAlgorithm}, f::MPIFile, params::CommonPreProcessingParameters{NoBackgroundCorrectionParameters})
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))), ignore = [:neglectBGFrames, :bgParams])
  result = getMeasurementsFD(f; bgCorrection = false, kwargs...)
  return result
end
function process(t::Type{<:AbstractMPIRecoAlgorithm}, f::MPIFile, params::CommonPreProcessingParameters{InternalBackgroundCorrectionParameters})
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))), ignore = [:neglectBGFrames, :bgParams])
  result = getMeasurementsFD(f; bgCorrection = true, interpolateBG = params.bgCorrection.interpolateBG, kwargs...)
  return result
end
function process(t::Type{<:AbstractMPIRecoAlgorithm}, f::MPIFile, params::CommonPreProcessingParameters{<:ExternalBackgroundCorrection})
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))), ignore = [:neglectBGFrames, :bgParams])
  result = getMeasurementsFD(f, bgCorrection = false, kwargs...)
  bgParams = fromKwargs(ExternalPreProcessedBackgroundCorrectionParameters; kwargs..., bgParams = params.bgParams)
  return process(t, result, bgParams)
end