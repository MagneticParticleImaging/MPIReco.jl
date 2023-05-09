export CommonPreProcessingParameters
Base.@kwdef struct CommonPreProcessingParameters{T<:AbstractBackgroundCorrectionParameters} <: AbstractPreProcessingParameters
  bgCorrection::T = NoBackgroundCorrection()
  neglectBGFrames::Bool = true
  numAverages::Int64 = 1
  numPeriodAverages::Int64 = 1
  numPeriodGrouping::Int64 = 1
  frames::Union{Nothing, UnitRange{Int64}} = nothing
  spectralLeakageCorrection::Bool = false
  loadasreal::Bool = false
end

function RecoUtils.process(t::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::CommonPreProcessingParameters{NoBackgroundCorrection})
  kwargs = toKwargs(params)
  if isnothing(params.frames)
    kwargs[:frames] = params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))
  end
  delete!(kwargs, :neglectBGFrames)
  result = getMeasurementsFD(f; bgCorrection = false, kwargs...)
  return result
end
function RecoUtils.process(t::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::CommonPreProcessingParameters{InternalBackgroundCorrection})
  kwargs = toKwargs(params)
  if isnothing(params.frames)
    kwargs[:frames] = params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))
  end
  delete!(kwargs, :neglectBGFrames)
  result = getMeasurementsFD(f; bgCorrection = true, interpolateBG = params.bgCorrection.interpolateBG, kwargs...)
  return result
end
function RecoUtils.process(t::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::CommonPreProcessingParameters{<:ExternalBackgroundCorrection})
  kwargs = toKwargs(params)
  if isnothing(params.frames)
    kwargs[:frames] = params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))
  end
  delete!(kwargs, :neglectBGFrames)
  result = getMeasurementsFD(f, bgCorrection = false, kwargs...)
  bgParams = fromKwargs(ExternalPreProcessedBackgroundCorrectionParameters; kwargs..., bgParams = params.bgCorrection)
  return process(t, result, bgParams)
end
RecoUtils.process(algo::AbstractMPIReconstructionAlgorithm, f::MPIFile, params::CommonPreProcessingParameters{<:AbstractBackgroundCorrectionParameters}) = process(typeof(algo), f, params)