Base.@kwdef struct CommonPreProcessingParameters{T<:AbstractBackgroundCorrectionParameters} <: AbstractPreProcessingParameter
  bgCorrection::T = NoBackgroundCorrection()
  numAverages::Int64 = 1
  numPeriodAverages::Int64 = 1
  numPeriodGrouping::Int64 = 1
  frames::Union{Nothing, UnitRange{Int64}} = nothing
  spectralLeakageCorrection::Bool = false
  loadasreal::Bool = false
end

# Not yet happy with field replication, but nested structs with parametric types is also quite ugly
Base.@kwdef struct FrequencyFilteredPreProcessingParamters{T<:AbstractBackgroundCorrectionParameters} <: AbstractPreProcessingParameter
  freqs::Vector{Int64}
  bgCorrection::T = NoBackgroundCorrection()
  numAverages::Int64 = 1
  numPeriodAverages::Int64 = 1
  numPeriodGrouping::Int64 = 1
  frames::Union{Nothing, UnitRange{Int64}} = nothing
  spectralLeakageCorrection::Bool = false
  loadasreal::Bool = false
end

function process(::Type{AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParamters{NoBackgroundCorrection})
  kwargs = toKwargs(params)
  isnothing(params.frames) && kwargs[:frames] = 1:acqNumFrames(f)
  result = getMeasurementFD(f, bgCorrection = false, frequencies = params.freqs, kwargs...)
  return result
end
function process(::Type{AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParamters{InternalBackgroundCorrection})
  kwargs = toKwargs(params)
  isnothing(params.frames) && kwargs[:frames] = 1:acqNumFrames(f)
  result = getMeasurementFD(f, bgCorrection = true, interpolateBG = params.bgCorrection.interpolateBG, frequencies = params.freqs, kwargs...)
  return result
end
function process(::Type{AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParamters{SimpleExternalBackgroundCorrection})
  kwargs = toKwargs(params)
  isnothing(params.frames) && kwargs[:frames] = 1:acqNumFrames(f)
  result = getMeasurementFD(f, bgCorrection = false, frequencies = params.freqs, kwargs...)
  bgParams = params.bgCorrection
  kwargs[:frames] = bgParams.bgFrames
  empty = getMeasurementsFD(bgParams.emptyMeas, false, bgCorrection = false, numAverages=length(bgFrames), frequencies = params.freqs, kwargs...)
  return result .-empty
end
function process(::Type{AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParamters{LinearInterpolatedExternalBackgroundCorrection})
  kwargs = toKwargs(params)
  isnothing(params.frames) && kwargs[:frames] = 1:acqNumFrames(f)
  result = getMeasurementFD(f, bgCorrection = false, frequencies = params.freqs, kwargs...)
  bgParams = params.bgCorrection
  kwargs[:frames] = bgParams.bgFrames
  empty = getMeasurementsFD(bgParams.emptyMeas, false, bgCorrection = false, numAverages=length(bgParams.bgFrames), frequencies = params.freqs, kwargs...)
  kwargs[:frames] = bgParams.bgFramesPost
  emptyPost =getMeasurementsFD(bgParams.emptyMeas, false, bgCorrection = false, numAverages=length(bgParams.bgFramesPost), frequencies = params.freqs, kwargs...)
  for l=1:size(result, 4)
    alpha = (l - mean(bgParams.bgFrames)) / (mean(bgParams.bgFramesPost) - mean(bgParams.bgFrames))
    result[:,:,l] .-=  (1-alpha).*empty[:,:,1] .+ alpha.*emptyPost[:,:,1]
  end
  return result
end
