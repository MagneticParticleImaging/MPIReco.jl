export SinglePatchReconstructionAlgorithm, SinglePatchReconstructionParameter, SinglePatchParameters

abstract type AbstractSinglePatchReconstructionAlgorithm <: AbstractMPIRecoAlgorithm end
abstract type AbstractSinglePatchReconstructionParameters <: AbstractMPIReconstructionParameters end
abstract type AbstractSinglePatchAlgorithmParameters <: AbstractMPIRecoParameters end

Base.@kwdef mutable struct SinglePatchParameters{PR<:AbstractMPIPreProcessingParameters,
     R<:AbstractSinglePatchReconstructionParameters, PT<:AbstractMPIPostProcessingParameters} <: AbstractSinglePatchAlgorithmParameters
  pre::Union{PR, AbstractUtilityReconstructionParameters{PR}}
  reco::R
  post::PT = NoPostProcessing() 
end
  
function (params::SinglePatchParameters)(algo::T, data::MPIFile, frequencies::Union{Vector{CartesianIndex{2}}, Nothing} = nothing) where {T<:AbstractSinglePatchReconstructionAlgorithm}
  result = params.pre(algo, data, frequencies)
  result = params.reco(algo, result)
  result = params.post(algo, result)
  result = finalizeResult(algo, result, data)
  return result
end
function (params::SinglePatchParameters)(algo::T, data::AbstractArray, args...) where {T<:AbstractSinglePatchReconstructionAlgorithm}
  throw(ArgumentError("SinglePatchAlgorithms are not defined for the given arguments, expected <: MPIFile, found $(typeof(data))"))
end

function finalizeResult(algo::AbstractSinglePatchReconstructionAlgorithm, result, data::MPIFile)
  pixspacing, offset = calcSpacingAndOffset(algo.sf, data, algo.grid)
  dt = acqNumAverages(data)*dfCycle(data)*numAverages(algo.params.pre)*1u"s"
  im = makeAxisArray(result, pixspacing, offset, dt)
  return ImageMeta(im, generateHeaderDict(algo.sf, data))
end

### Raw Data
#consistenceCheck(algo.sf::MPIFile, u::AbstractArray{ComplexF})


include("SinglePatchAlgorithm.jl")
include("SinglePatchTwoStepReconstructionAlgorithm.jl")
include("SinglePatchBGEstimationAlgorithm.jl")
#include("SinglePatchTemporalRegularizationAlgorithm.jl")

### Multi-Threading
#=
consistenceCheck(sf::MPIFile, threaded::MultiThreadedInput) = consistenceCheck(sf, threaded.inputs[1]) 
finalizeResult(algo::AbstractSinglePatchReconstructionAlgorithm, result, threadedInput::MultiThreadedInput) = finalizeResult(algo, result, threadedInput.inputs[1])

function (params::SinglePatchParameters)(algo::T, threadedInput::MultiThreadedInput, frequencies::Union{Vector{CartesianIndex{2}}, Nothing} = nothing) where {T<:AbstractSinglePatchReconstructionAlgorithm}
  result = params.pre(algo, threadedInput, frequencies)
  result = params.reco(algo, MultiThreadedInput(threadedInput.scheduler, (result,)))
  result = params.post(algo, MultiThreadedInput(threadedInput.scheduler, (result,)))
  return result
end

function (params::P)(algo::T, threadInput::MultiThreadedInput, frequencies::Vector{CartesianIndex{2}}) where {T<:Union{SinglePatchReconstructionAlgorithm, SinglePatchBGEstimationAlgorithm}, P <: AbstractMPIPreProcessingParameters}
  scheduler = threadInput.scheduler
  data = threadInput.inputs

  # Compute how many blocks of frames are averaged into one reco frame
  paramFrames = params.frames
  numFrames = length(paramFrames)
  averages = params.numAverages
  start = collect(1:averages:numFrames)
  stop = map(x->min(x+averages-1, numFrames), start)
  blocks = zip(start, stop)

  # Distribute frames onto nthreads tasks
  numThreads = AbstractImageReconstruction.nthreads(scheduler)
  threadFrames = Vector{Vector{Int64}}()
  for block in Iterators.partition(blocks, div(length(blocks), numThreads))
    tmp = map(x-> collect(paramFrames[x[1]:x[2]]), block)
    push!(threadFrames, vcat(tmp...))
  end

  for frames in threadFrames
    pre = fromKwargs(P; toKwargs(params, flatten = DataType[])..., frames = frames)
    put!(scheduler, algo, pre, data..., frequencies)
  end
  result = nothing
  for frames in threadFrames
    if isnothing(result)
      result = take!(scheduler)
    else
      result = cat(result, take!(scheduler); dims = 3)
    end
  end
  return result
end

function (params::LeastSquaresParameters)(algo::T, threadInput::MultiThreadedInput) where {T<:Union{SinglePatchReconstructionAlgorithm, SinglePatchBGEstimationAlgorithm}}
  scheduler = threadInput.scheduler
  data = threadInput.inputs
  u = threadInput.inputs[1]
  L = size(u)[end]

  # Distribute frames onto nthreads tasks
  numThreads = AbstractImageReconstruction.nthreads(scheduler)
  threadFrames = Vector{UnitRange{Int64}}()
  for block in Iterators.partition(1:L, ceil(Int64, L/numThreads))
    push!(threadFrames, block)
  end

  for frames in threadFrames
    put!(scheduler, algo, params, reshape(u, :, L)[:,frames])
  end
  result = nothing
  for frames in threadFrames
    if isnothing(result)
      result = take!(scheduler)
    else
      result = cat(result, take!(scheduler); dims = 2)
    end
  end
  return result
end

(params::NoPostProcessing)(algo::T, threadInput::MultiThreadedInput) where {T<:Union{SinglePatchReconstructionAlgorithm, SinglePatchBGEstimationAlgorithm}} = params(algo, threadInput.inputs...)
=#