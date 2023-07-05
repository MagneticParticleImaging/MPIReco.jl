export AbstractFrequencyFilterParameter
abstract type AbstractFrequencyFilterParameter <: AbstractMPIRecoParameters end
# Could possible also be nested
export SNRThresholdFrequencyFilterParameter
Base.@kwdef struct SNRThresholdFrequencyFilterParameter <: AbstractFrequencyFilterParameter
  minFreq::Float64 = -1.0
  maxFreq::Float64 = -1.0
  recChannels::UnitRange{Int64} = 1:1
  SNRThresh::Float64=-1.0
  sortBySNR::Bool = false
  numPeriodAverages::Int64 = 1
  numPeriodGrouping::Int64 = 1
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
  maxFreq::Float64 = 0.0
  recChannels::UnitRange{Int64} = 1:1
  numUsedFreqs::Int64=1
  sortBySNR::Bool = false
  numPeriodAverages::Int64 = 1
  numPeriodGrouping::Int64 = 1
end
RecoUtils.process(::Type{<:AbstractMPIReconstructionAlgorithm}, file::MPIFile, params::AbstractFrequencyFilterParameter) = filterFrequencies(file; toKwargs(params)...)
