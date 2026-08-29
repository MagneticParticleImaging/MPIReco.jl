export AbstractFrequencyFilterParameter
abstract type AbstractFrequencyFilterParameter <: AbstractMPIRecoParameters end

export NoFrequencyFilterParameter
# TODO This requires numPeriodAverages and numPeriodGrouping and should use a freq. loading function from MPIFiles
struct NoFrequencyFilterParameter <: AbstractFrequencyFilterParameter end

function (params::NoFrequencyFilterParameter)(::Type{<:AbstractMPIRecoAlgorithm}, file::MPIFile)
  return collect(vec(CartesianIndices((rxNumFrequencies(file), rxNumChannels(file)))))
end

export DirectSelectionFrequencyFilterParameters
struct DirectSelectionFrequencyFilterParameters{T <: Union{Integer, CartesianIndex{2}}, FIT <: AbstractVector{T}} <: AbstractFrequencyFilterParameter
  freqIndices::Union{Nothing, FIT}
  function DirectSelectionFrequencyFilterParameters(;freqIndices::FIT = nothing) where FIT
    if isnothing(freqIndices)
      el = CartesianIndex{2}
      v = Vector{el}
    else
      el = eltype(freqIndices)
      v = typeof(freqIndices)
    end
    new{el, v}(freqIndices)
  end
end
function (params::DirectSelectionFrequencyFilterParameters)(::Type{<:AbstractMPIRecoAlgorithm}, file::MPIFile) where T <: Integer
  nFreq = params.freqIndices
  nReceivers = rxNumChannels(file)
  return vec([CartesianIndex{2}(i, j) for i in nFreq, j in nReceivers])
end
function (params::DirectSelectionFrequencyFilterParameters{T})(::Type{<:AbstractMPIRecoAlgorithm}, file::MPIFile) where T <: CartesianIndex{2}
  return params.freqIndices
end

# Could possible also be nested
export SNRThresholdFrequencyFilterParameter
@parameter struct SNRThresholdFrequencyFilterParameter <: AbstractFrequencyFilterParameter
  minFreq::Float64 = 0.0
  maxFreq::Union{Float64, Nothing} = nothing
  recChannels::Union{Vector{Int64}, UnitRange{Int64}, Nothing} = nothing
  SNRThresh::Float64=-1.0
  #sortBySNR::Bool = false
  numPeriodAverages::Int64 = 1
  numPeriodGrouping::Int64 = 1
  maxMixingOrder::Int64 = -1
  numSidebandFreqs::Int64 = -1
end
export defaultParameterMaxFreq, defaultParameterMinFreq, defaultParameterRecChannels
defaultParameterMaxFreq(new::MPIFile) = rxBandwidth(new)
defaultParameterMaxFreq(new::Missing) = missing
defaultParameterMinFreq(new::MPIFile) = first(rxFrequencies(new))
defaultParameterMinFreq(new::Missing) = missing
defaultParameterRecChannels(new::MPIFile) = 1:rxNumChannels(new)
defaultParameterRecChannels(new::Missing) = missing
#function SNRThresholdFrequencyFilterParameter(;sf::MPIFile, maxFreq = nothing, recChannels = nothing)
#  if isnothing(maxFreq)
#    maxFreq = rxBandwidth(sf)
#  end
#  if isnothing(recChannels)
#    recChannels = 1:rxNumChannels(sf)
#  end
#  return SNRThresholdFrequencyFilterParameter(;maxFreq = maxFreq, recChannels = recChannels)
#end
export FreqNumThresholdFrequencyFilterParameter
@parameter struct FreqNumThresholdFrequencyFilterParameter <: AbstractFrequencyFilterParameter
  minFreq::Float64 = 0.0
  maxFreq::Union{Float64, Nothing} = nothing
  recChannels::Union{UnitRange{Int64}, Nothing} = nothing
  numUsedFreqs::Int64=1
  #sortBySNR::Bool = false
  numPeriodAverages::Int64 = 1
  numPeriodGrouping::Int64 = 1
  maxMixingOrder::Int64 = -1
  numSidebandFreqs::Int64 = -1
end

function (params::AbstractFrequencyFilterParameter)(::Type{<:AbstractMPIRecoAlgorithm}, file::MPIFile)
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:maxFreq => rxBandwidth(file), :recChannels => 1:rxNumChannels(file))) 
  filterFrequencies(file; kwargs...)
end

export CompositeFrequencyFilterParameters
@parameter struct CompositeFrequencyFilterParameters <: AbstractFrequencyFilterParameter
  filters::Vector{AbstractFrequencyFilterParameter}
end
function (params::CompositeFrequencyFilterParameters)(algoT::Type{<:AbstractMPIRecoAlgorithm}, file::MPIFile)
  return reduce(intersect, filter(!isnothing, map(p -> p(algoT, file), params.filters)))
end

#=
export NoiseLevelFrequencyFilterParameter
@parameter struct NoiseLevelFrequencyFilterParameter <: AbstractFrequencyFilterParameter
  noiseMeas::MPIFile
  levelFactor::Float64
  noiseFrames::Union{Vector{Int64}, UnitRange{Int64}} = measBGFrameIdx(noiseMeas)
  minFreq::Float64 = 0.0
  maxFreq::Union{Float64, Nothing} = nothing
  recChannels::Union{UnitRange{Int64}, Nothing} = nothing
  numPeriodAverages::Int64 = 1
  numPeriodGrouping::Int64 = 1
  maxMixingOrder::Int64 = -1
  numSidebandFreqs::Int64 = -1
end
function (params::NoiseLevelFrequencyFilterParameter)(::Type{<:AbstractMPIRecoAlgorithm}, file::MPIFile)
  noiseLevel = getNoiseLevel(params.noiseMeas, params.frames, params.recChannels)
  kwargs = toKwargs(params, ignore = [:noiseMeas, :frames], default = Dict{Symbol, Any}(:maxFreq => rxBandwidth(file), :recChannels => 1:rxNumChannels(file)))
  return filterFrequencies(file; kwargs..., SNRThresh = params.levelFactor * noiseLevel)
end=#