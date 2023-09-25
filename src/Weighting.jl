export getWeights, WeightingType, setFreqToZero

export AbstractWeightingParameters
abstract type AbstractWeightingParameters <: AbstractMPIRecoParameters end

export NoWeightingParamters
struct NoWeightingParamters <: AbstractWeightingParameters end
process(::Type{<:AbstractMPIRecoAlgorithm}, params::NoWeightingParamters, data) = nothing

export ChannelWeightingParameters
Base.@kwdef struct ChannelWeightingParameters <: AbstractWeightingParameters
  channelWeights::Vector{Float64} = [1.0, 1.0, 1.0]
end
process(::Type{<:AbstractMPIRecoAlgorithm}, params::ChannelWeightingParameters, data::Vector{CartesianIndex{2}}) = map(x-> params.channelWeights[x[2]], data)

#=
baremodule WeightingType
  None = 0
  Norm = 1
  MixingFactors = 2
  BGVariance = 3
  VarMeanThresh = 4
  Channel = 5
end

function getWeights(weightType, freq, S; weightingLimit=0.0, emptyMeas = nothing,
                 bgFrames=1:10, bMeas = nothing, fgFrames = 1:10, bSF=nothing,
                 channelWeights=[1.0,1.0,1.0])

  if weightType == WeightingType.None
    return nothing
  elseif weightType == WeightingType.Norm
    reciprocalWeights = rowEnergy(S)
  elseif weightType == WeightingType.MixingFactors
    error("This weighting mode has to be implemented")
  elseif weightType == WeightingType.Channel
    return map(f -> channelWeights[f[2]], freq)
  elseif weightType == WeightingType.BGVariance
    if emptyMeas == nothing
      stdDevU = sqrt.(vec(getBV(bSF)))
    else
      uEmpty = getMeasurementsFT(emptyMeas,frames=bgFrames)
      stdDevU = sqrt(abs(var(uEmpty,3 )))
    end
    reciprocalWeights = stdDevU[freq]
  elseif weightType == WeightingType.VarMeanThresh

    nFreq = numFreq(bMeas)
    nReceivers = numReceivers(bMeas)
    freqMask = zeros(nFreq, nReceivers)

    u = getMeasurementsFT(bMeas,frames=fgFrames, spectralCleaning=true)
    uEmpty = getMeasurementsFT(emptyMeas,frames=bgFrames, spectralCleaning=true)

    stdDevU = sqrt(abs(var(u,3 )))
    stdDevUEmpty = sqrt(abs(var(uEmpty,3 )))
    meanU = abs(mean(u,dims=3).-mean(uEmpty,dims=3))
    meanUEmpty = abs(mean(uEmpty, dims=3))

    for k=1:nFreq
      for r=1:nReceivers
        if stdDevU[k,r]/meanU[k,r] < 0.5
          freqMask[k,r] = 1
        else
          freqMask[k,r] = 1000
        end
      end
    end

    reciprocalWeights = vec(freqMask)[freq]
  end

  if weightingLimit>0

     m = maximum(reciprocalWeights)

     # The idea here is to only normalize those rows which have enough energy
     reciprocalWeights[ reciprocalWeights .< m*weightingLimit ] = m*weightingLimit
  end

  weights = copy(reciprocalWeights)

  for l=1:length(weights)
    weights[l] = 1 / reciprocalWeights[l].^2
  end

  return weights
end


function setNoiseFreqToZero(uMeas, freq, noiseFreqThresh; emptyMeas = nothing, bgFrames=1:10, bMeas = nothing, fgFrames = 1:10)
  @debug "Setting noise frequencies to zero"

  nFreq = numFreq(bMeas)
  nReceivers = numReceivers(bMeas)
  freqMask = zeros(nFreq, nReceivers)

  u = getMeasurementsFT(bMeas,frames=fgFrames)
  uEmpty = getMeasurementsFT(emptyMeas,frames=bgFrames)

  stdDevU = sqrt(abs(var(u,3 )))
  meanU = abs(mean(u,dims=3).-mean(uEmpty,dims=3))
  meanUEmpty = abs(mean(uEmpty,dims=3))

  for k=1:nFreq
    for r=1:nReceivers
      if stdDevU[k,r]/meanU[k,r] < noiseFreqThresh
        freqMask[k,r] = 1
      else
        freqMask[k,r] = 0
      end
    end
  end

  uMeas[:,:] .*=  vec(freqMask)[freq]
end
=#