export AbstractFrequencyFilterParameter
abstract type AbstractFrequencyFilterParameter <: AbstractMPIRecoParameters end
# Could possible also be nested
export SNRThresholdFrequencyFilterParameter
Base.@kwdef struct SNRThresholdFrequencyFilterParameter <: AbstractFrequencyFilterParameter
  minFreq::Float64 = 0.0
  maxFreq::Union{Float64, Nothing} = nothing
  recChannels::Union{UnitRange{Int64}, Nothing} = nothing
  SNRThresh::Float64=-1.0
  #sortBySNR::Bool = false
  numPeriodAverages::Int64 = 1
  numPeriodGrouping::Int64 = 1
  maxMixingOrder::Int64 = -1
  numSidebandFreqs::Int64 = -1
end
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
Base.@kwdef struct FreqNumThresholdFrequencyFilterParameter <: AbstractFrequencyFilterParameter
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
function process(::Type{<:AbstractMPIRecoAlgorithm}, params::AbstractFrequencyFilterParameter, file::MPIFile)
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:maxFreq => rxBandwidth(file), :recChannels => 1:rxNumChannels(file))) 
  filterFrequencies(file; kwargs...)
end
