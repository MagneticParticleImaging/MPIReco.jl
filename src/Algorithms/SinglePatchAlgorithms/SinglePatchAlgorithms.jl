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
  
function process(algo::T, params::SinglePatchParameters, data::MPIFile, frequencies::Union{Vector{CartesianIndex{2}}, Nothing} = nothing) where {T<:AbstractSinglePatchReconstructionAlgorithm}
  result = process(algo, params.pre, data, frequencies)
  result = process(algo, params.reco, result)
  result = process(algo, params.post, result)
  return result
end

function AbstractImageReconstruction.put!(algo::AbstractSinglePatchReconstructionAlgorithm, data)
  lock(algo) do
    #consistenceCheck(algo.sf, data)
    
    result = process(algo, algo.params, data, algo.freqs)

    result = finalizeResult(algo, result, data)

    Base.put!(algo.output, result)
  end
end

function finalizeResult(algo::AbstractSinglePatchReconstructionAlgorithm, result, data::MPIFile)
  pixspacing = (spacing(algo.grid) ./ acqGradient(data)[1] .* acqGradient(algo.sf)[1])*1000u"mm"
  offset = (ffPos(data) .- 0.5 .* calibFov(algo.sf))*1000u"mm" .+ 0.5 .* pixspacing
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

function process(algo::T, params::SinglePatchParameters, threadedInput::MultiThreadedInput, frequencies::Union{Vector{CartesianIndex{2}}, Nothing} = nothing) where {T<:AbstractSinglePatchReconstructionAlgorithm}
  result = process(algo, params.pre, threadedInput, frequencies)
  result = process(algo, params.reco, MultiThreadedInput(threadedInput.scheduler, (result,)))
  result = process(algo, params.post, MultiThreadedInput(threadedInput.scheduler, (result,)))
  return result
end

function process(algo::T, params::P, threadInput::MultiThreadedInput, frequencies::Vector{CartesianIndex{2}}) where {T<:Union{SinglePatchReconstructionAlgorithm, SinglePatchBGEstimationAlgorithm}, P <: AbstractMPIPreProcessingParameters}
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

function process(algo::T, params::LeastSquaresParameters, threadInput::MultiThreadedInput) where {T<:Union{SinglePatchReconstructionAlgorithm, SinglePatchBGEstimationAlgorithm}}
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

process(algo::T, params::NoPostProcessing, threadInput::MultiThreadedInput) where {T<:Union{SinglePatchReconstructionAlgorithm, SinglePatchBGEstimationAlgorithm}} = process(algo, params, threadInput.inputs...)
=#