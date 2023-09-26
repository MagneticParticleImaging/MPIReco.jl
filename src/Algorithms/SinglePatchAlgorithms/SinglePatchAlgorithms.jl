export SinglePatchReconstructionAlgorithm, SinglePatchReconstructionParameter, SinglePatchParameters

abstract type AbstractSinglePatchReconstructionAlgorithm <: AbstractMPIRecoAlgorithm end
abstract type AbstractSinglePatchReconstructionParameters <: AbstractMPIReconstructionParameters end
abstract type AbstractSinglePatchAlgorithmParameters <: AbstractMPIRecoParameters end

Base.@kwdef mutable struct SinglePatchParameters{PR<:AbstractMPIPreProcessingParameters,
     R<:AbstractSinglePatchReconstructionParameters, PT<:AbstractMPIPostProcessingParameters} <: AbstractSinglePatchAlgorithmParameters
  pre::PR
  reco::R
  post::PT = NoPostProcessing() 
end
  
function process(algo::T, params::SinglePatchParameters, data::MPIFile, frequencies::Union{Vector{CartesianIndex{2}}, Nothing} = nothing) where {T<:AbstractSinglePatchReconstructionAlgorithm}
  result = process(algo, params.pre, data, frequencies)
  result = process(algo, params.reco, result)
  result = process(algo, params.post, result)
  return result
end

function AbstractImageReconstruction.put!(algo::AbstractSinglePatchReconstructionAlgorithm, data)
  consistenceCheck(algo.sf, data)
  
  result = process(algo, algo.params, data, algo.freqs)

  result = finalizeResult(algo, result, data)

  Base.put!(algo.output, result)
end

function finalizeResult(algo::AbstractSinglePatchReconstructionAlgorithm, result, data::MPIFile)
  pixspacing = (spacing(algo.grid) ./ acqGradient(data)[1] .* acqGradient(algo.sf)[1])*1000u"mm"
  offset = (ffPos(data) .- 0.5 .* calibFov(algo.sf))*1000u"mm" .+ 0.5 .* pixspacing
  dt = acqNumAverages(data)*dfCycle(data)*algo.params.pre.numAverages*1u"s"
  im = makeAxisArray(result, pixspacing, offset, dt)
  return ImageMeta(im, generateHeaderDict(algo.sf, data))
end


include("SinglePatchAlgorithm.jl")
include("SinglePatchTwoStepReconstructionAlgorithm.jl")
include("SinglePatchBGEstimationAlgorithm.jl")
#include("SinglePatchTemporalRegularizationAlgorithm.jl")

### Multi-Threading
consistenceCheck(sf::MPIFile, threaded::MultiThreadedInput) = consistenceCheck(sf, threaded.inputs[1]) 
finalizeResult(algo::AbstractSinglePatchReconstructionAlgorithm, result, threadedInput::MultiThreadedInput) = finalizeResult(algo, result, threadedInput.inputs[1])
function process(algo::T, params::SinglePatchParameters, threadInput::MultiThreadedInput, frequencies::Vector{CartesianIndex{2}}) where {T<:Union{SinglePatchReconstructionAlgorithm, SinglePatchTwoStepReconstructionAlgorithm, SinglePatchBGEstimationAlgorithm}}
  scheduler = threadInput.scheduler
  data = threadInput.inputs
  frames = collect(params.pre.frames) 
  # TODO consider averaging and periodgrouping
  for frame in frames
    pre = fromKwargs(CommonPreProcessingParameters; toKwargs(params.pre)..., frames = [frame])
    p = SinglePatchParameters(pre = pre, reco = params.reco, post = params.post)
    put!(scheduler, algo, p, data..., frequencies)
  end
  result = nothing
  meta = nothing
  for frame in frames
    if isnothing(result)
      result = take!(scheduler)
      meta = result
    else
      result = cat(result, take!(scheduler); dims = 5)
    end
  end
  return result
end