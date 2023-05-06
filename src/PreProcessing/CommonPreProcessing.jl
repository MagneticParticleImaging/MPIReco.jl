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

function RecoUtils.process(::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::CommonPreProcessingParameters{NoBackgroundCorrection})
  kwargs = toKwargs(params)
  if isnothing(params.frames)
    kwargs[:frames] = params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))
  end
  delete!(kwargs, :neglectBGFrames)
  result = getMeasurementsFD(f; bgCorrection = false, kwargs...)
  return result
end
function RecoUtils.process(::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::CommonPreProcessingParameters{InternalBackgroundCorrection})
  kwargs = toKwargs(params)
  if isnothing(params.frames)
    kwargs[:frames] = params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))
  end
  delete!(kwargs, :neglectBGFrames)
  result = getMeasurementsFD(f; bgCorrection = true, interpolateBG = params.bgCorrection.interpolateBG, kwargs...)
  return result
end
function RecoUtils.process(::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::CommonPreProcessingParameters{SimpleExternalBackgroundCorrection})
  kwargs = toKwargs(params)
  if isnothing(params.frames)
    kwargs[:frames] = params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))
  end
  delete!(kwargs, :neglectBGFrames)
  result = getMeasurementsFD(f, bgCorrection = false, kwargs...)
  bgParams = params.bgCorrection
  kwargs[:frames] = bgParams.bgFrames
  empty = getMeasurementsFD(bgParams.emptyMeas, false, bgCorrection = false, numAverages=length(bgFrames), kwargs...)
  return result .-empty
end
function RecoUtils.process(::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::CommonPreProcessingParameters{LinearInterpolatedExternalBackgroundCorrection})
  kwargs = toKwargs(params)
  if isnothing(params.frames)
    kwargs[:frames] = params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))
  end
  delete!(kwargs, :neglectBGFrames)
  result = getMeasurementsFD(f, bgCorrection = false, kwargs...)
  bgParams = params.bgCorrection
  kwargs[:frames] = bgParams.bgFrames
  empty = getMeasurementsFD(bgParams.emptyMeas, false, bgCorrection = false, numAverages=length(bgParams.bgFrames), kwargs...)
  kwargs[:frames] = bgParams.bgFramesPost
  emptyPost = getMeasurementsFD(bgParams.emptyMeas, false, bgCorrection = false, numAverages=length(bgParams.bgFramesPost), kwargs...)
  for l=1:size(result, 4)
    alpha = (l - mean(bgParams.bgFrames)) / (mean(bgParams.bgFramesPost) - mean(bgParams.bgFrames))
    result[:,:,l] .-=  (1-alpha).*empty[:,:,1] .+ alpha.*emptyPost[:,:,1]
  end
  return result
end
RecoUtils.process(algo::AbstractMPIReconstructionAlgorithm, f::MPIFile, params::CommonPreProcessingParameters{<:AbstractBackgroundCorrectionParameters}) = process(typeof(algo), f, params)