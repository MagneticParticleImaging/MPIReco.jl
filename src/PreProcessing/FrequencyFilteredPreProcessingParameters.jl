# Not yet happy with field replication, but nested structs with parametric types is also quite ugly
export FrequencyFilteredPreProcessingParameters
Base.@kwdef struct FrequencyFilteredPreProcessingParameters{T<:AbstractBackgroundCorrectionParameters} <: AbstractPreProcessingParameters
  freqs::Vector{Int64}
  neglectBGFrames::Bool = true
  bgCorrection::T = NoBackgroundCorrection()
  numAverages::Int64 = 1
  numPeriodAverages::Int64 = 1
  numPeriodGrouping::Int64 = 1
  frames::Union{Nothing, UnitRange{Int64}} = nothing
  spectralLeakageCorrection::Bool = false
  loadasreal::Bool = false
end

function RecoUtils.process(::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParameters{NoBackgroundCorrection})
  kwargs = toKwargs(params)
  if isnothing(params.frames)
    kwargs[:frames] = params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))
  end
  delete!(kwargs, :neglectBGFrames)
  result = getMeasurementsFD(f, bgCorrection = false, frequencies = params.freqs, kwargs...)
  return result
end
function RecoUtils.process(::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParameters{InternalBackgroundCorrection})
  kwargs = toKwargs(params)
  if isnothing(params.frames)
    kwargs[:frames] = params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))
  end
  delete!(kwargs, :neglectBGFrames)
  result = getMeasurementsFD(f, bgCorrection = true, interpolateBG = params.bgCorrection.interpolateBG, frequencies = params.freqs, kwargs...)
  return result
end
function RecoUtils.process(::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParameters{SimpleExternalBackgroundCorrection})
  kwargs = toKwargs(params)
  if isnothing(params.frames)
    kwargs[:frames] = params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))
  end
  delete!(kwargs, :neglectBGFrames)
  result = getMeasurementFD(f, bgCorrection = false, frequencies = params.freqs, kwargs...)
  bgParams = params.bgCorrection
  kwargs[:frames] = bgParams.bgFrames
  empty = getMeasurementsFD(bgParams.emptyMeas, false, bgCorrection = false, numAverages=length(bgFrames), frequencies = params.freqs, kwargs...)
  return result .-empty
end
function RecoUtils.process(::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParameters{LinearInterpolatedExternalBackgroundCorrection})
  kwargs = toKwargs(params)
  if isnothing(params.frames)
    kwargs[:frames] = params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))
  end
  delete!(kwargs, :neglectBGFrames)
  result = getMeasurementsFD(f, bgCorrection = false, frequencies = params.freqs, kwargs...)
  bgParams = params.bgCorrection
  kwargs[:frames] = bgParams.bgFrames
  empty = getMeasurementsFD(bgParams.emptyMeas, false, bgCorrection = false, numAverages=length(bgParams.bgFrames), frequencies = params.freqs, kwargs...)
  kwargs[:frames] = bgParams.bgFramesPost
  emptyPost = getMeasurementsFD(bgParams.emptyMeas, false, bgCorrection = false, numAverages=length(bgParams.bgFramesPost), frequencies = params.freqs, kwargs...)
  for l=1:size(result, 4)
    alpha = (l - mean(bgParams.bgFrames)) / (mean(bgParams.bgFramesPost) - mean(bgParams.bgFrames))
    result[:,:,l] .-=  (1-alpha).*empty[:,:,1] .+ alpha.*emptyPost[:,:,1]
  end
  return result
end
RecoUtils.process(algo::AbstractMPIReconstructionAlgorithm, f::MPIFile, params::FrequencyFilteredPreProcessingParameters{<:AbstractBackgroundCorrectionParameters}) = process(typeof(algo), f, params)