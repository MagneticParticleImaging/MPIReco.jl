export MultiPatchReconstructionAlgorithm, MultiPatchReconstructionParameter
Base.@kwdef struct MultiPatchReconstructionParameter{arrT <: AbstractArray,F<:AbstractFrequencyFilterParameter,O<:AbstractMultiPatchOperatorParameter, S<:AbstractSolverParameters, FF<:AbstractFocusFieldPositions, FFSF<:AbstractFocusFieldPositions, R <: AbstractRegularization, W<:AbstractWeightingParameters} <: AbstractMultiPatchReconstructionParameters
  arrayType::Type{arrT} = Array
  # File
  sf::MultiMPIFile
  freqFilter::F
  opParams::Union{O, AbstractUtilityReconstructionParameters{O}}
  ffPos::FF = DefaultFocusFieldPositions()
  ffPosSF::FFSF = DefaultFocusFieldPositions()
  solverParams::S
  reg::Vector{R} = AbstractRegularization[]
  weightingParams::Union{W, AbstractUtilityReconstructionParameters{W}} = NoWeightingParameters()
end

Base.@kwdef mutable struct MultiPatchReconstructionAlgorithm{P, arrT <: AbstractArray, vecT <: AbstractArray} <: AbstractMultiPatchReconstructionAlgorithm where {P<:AbstractMultiPatchAlgorithmParameters}
  params::P
  # Could also do reconstruction progress meter here
  opParams::Union{AbstractMultiPatchOperatorParameter, AbstractUtilityReconstructionParameters{<:AbstractMultiPatchOperatorParameter},Nothing} = nothing
  sf::MultiMPIFile
  weights::Union{Nothing, vecT} = nothing
  arrayType::Type{arrT}
  ffOp::Union{Nothing, AbstractMultiPatchOperator}
  ffPos::Union{Nothing,AbstractArray}
  ffPosSF::Union{Nothing,AbstractArray}
  freqs::Vector{CartesianIndex{2}}
  output::Channel{Any}
end

function MultiPatchReconstruction(params::MultiPatchParameters{<:AbstractMPIPreProcessingParameters,R,PT}) where {R<:AbstractMultiPatchReconstructionParameters,PT<:AbstractMPIPostProcessingParameters}
  return MultiPatchReconstructionAlgorithm(params)
end
function MultiPatchReconstructionAlgorithm(params::MultiPatchParameters{<:AbstractMPIPreProcessingParameters,<:MultiPatchReconstructionParameter,<:AbstractMPIPostProcessingParameters})
  reco = params.reco
  freqs = process(MultiPatchReconstructionAlgorithm, reco.freqFilter, reco.sf)

  # Prepare operator construction
  ffPos_ = positions(reco.ffPos)
  ffPosSF = positions(reco.ffPosSF)
  if isnothing(ffPosSF)
    L = length(ffPos(reco.sf[1]))
    ffPosSF = [vec(ffPos(SF))[l] for l=1:L, SF in reco.sf]
  end

  return MultiPatchReconstructionAlgorithm{typeof(params), reco.arrayType, typeof(reco.arrayType{Float32}(undef, 0))}(params, reco.opParams, reco.sf, nothing, reco.arrayType, nothing, ffPos_, ffPosSF, freqs, Channel{Any}(Inf))
end
recoAlgorithmTypes(::Type{MultiPatchReconstruction}) = SystemMatrixBasedAlgorithm()
AbstractImageReconstruction.parameter(algo::MultiPatchReconstructionAlgorithm) = algo.origParam

Base.lock(algo::MultiPatchReconstructionAlgorithm) = lock(algo.output)
Base.unlock(algo::MultiPatchReconstructionAlgorithm) = unlock(algo.output)
Base.isready(algo::MultiPatchReconstructionAlgorithm) = isready(algo.output)
Base.wait(algo::MultiPatchReconstructionAlgorithm) = wait(algo.output)
AbstractImageReconstruction.take!(algo::MultiPatchReconstructionAlgorithm) = Base.take!(algo.output)

function AbstractImageReconstruction.put!(algo::MultiPatchReconstructionAlgorithm, data::MPIFile)
  lock(algo) do 
    #consistenceCheck(algo.sf, data)

    algo.ffOp, algo.weights = process(algo, algo.opParams, data, algo.freqs, algo.params.reco.weightingParams)
    
    result = process(algo, algo.params, data, algo.freqs)

  # Create Image (maybe image parameter as post params?)
  # TODO make more generic to apply to other pre/reco params as well (pre.numAverage main issue atm)
  pixspacing = (voxelSize(algo.sf) ./ sfGradient(data,3) .* sfGradient(algo.sf,3)) * 1000u"mm"
  offset = (fieldOfViewCenter(algo.ffOp.grid) .- 0.5.*fieldOfView(algo.ffOp.grid) .+ 0.5.*spacing(algo.ffOp.grid)) * 1000u"mm"
  dt = acqNumAverages(data) * dfCycle(data) * numAverages(algo.params.pre) * 1u"s"
  im = makeAxisArray(result, pixspacing, offset, dt)
  result = ImageMeta(im, generateHeaderDict(algo.sf, data))

    Base.put!(algo.output, result)
  end
end

function process(algo::MultiPatchReconstructionAlgorithm, params::Union{OP, ProcessResultCache{OP}}, f::MPIFile, frequencies::Vector{CartesianIndex{2}}, weightingParams) where OP <: AbstractMultiPatchOperatorParameter
  ffPos_ = ffPos(f)
  periodsSortedbyFFPos = unflattenOffsetFieldShift(ffPos_)
  idxFirstPeriod = getindex.(periodsSortedbyFFPos,1)
  ffPos_ = ffPos_[:,idxFirstPeriod]
  gradient = acqGradient(f)[:,:,1,idxFirstPeriod]

  if !isnothing(algo.ffPos)
    ffPos_[:] = algo.ffPos
  end
  
  result = process(typeof(algo), params, algo.sf, frequencies, gradient, ffPos_, algo.ffPosSF)
  # Kinda of hacky. MultiPatch parameters don't map nicely to the SinglePatch inspired pre, reco, post structure
  # Have to create weights before ffop is (potentially) moved to GPU, as GPU arrays don't have efficient hash implementations
  # Which makes this process expensive to cache
  weights = process(typeof(algo), weightingParams, frequencies, result, nothing, algo.arrayType)
  resultXPU = process(typeof(algo), params, result, algo.arrayType)
  return resultXPU, weights
end

function process(algo::MultiPatchReconstructionAlgorithm, params::Union{A, ProcessResultCache{<:A}}, f::MPIFile, args...) where A <: AbstractMPIPreProcessingParameters
  result = process(typeof(algo), params, f, args...)
  result = adapt(algo.arrayType, result)
  return result
end

function process(t::Type{<:MultiPatchReconstructionAlgorithm}, params::CommonPreProcessingParameters{NoBackgroundCorrectionParameters}, f::MPIFile, frequencies::Union{Vector{CartesianIndex{2}}, Nothing} = nothing)
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))), ignore = [:neglectBGFrames, :bgCorrection])
  result = getMeasurementsFD(f, bgCorrection = false; kwargs..., frequencies = frequencies)
  periodsSortedbyFFPos = unflattenOffsetFieldShift(ffPos(f))
  uTotal = Base.similar(result,size(result,1),length(periodsSortedbyFFPos),size(result,3))
  for k=1:length(periodsSortedbyFFPos)
      uTotal[:,k,:] = mean(result[:,periodsSortedbyFFPos[k],:], dims=2)
  end
  return uTotal
end

function process(t::Type{<:MultiPatchReconstructionAlgorithm}, params::CommonPreProcessingParameters{SimpleExternalBackgroundCorrectionParameters}, f::MPIFile, frequencies::Union{Vector{CartesianIndex{2}}, Nothing} = nothing)
  # Foreground, ignore BGCorrection to reuse preprocessing
  fgParams = CommonPreProcessingParameters(;toKwargs(params)..., bgParams = NoBackgroundCorrectionParameters())
  result = process(t, fgParams, f, frequencies)
  # Background
  kwargs = toKwargs(params, ignore = [:neglectBGFrames, :bgCorrection],
    default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))))
  bgParams = fromKwargs(ExternalPreProcessedBackgroundCorrectionParameters; kwargs..., bgParams = params.bgParams)
  return process(t, bgParams, result, frequencies)
end

function process(::Type{<:MultiPatchReconstructionAlgorithm}, params::ExternalPreProcessedBackgroundCorrectionParameters{SimpleExternalBackgroundCorrectionParameters}, data::Array, frequencies::Union{Vector{CartesianIndex{2}}, Nothing} = nothing)
  kwargs = toKwargs(params, overwrite = Dict{Symbol, Any}(:frames => params.bgParams.bgFrames), ignore = [:bgParams])
  # TODO migrate with hardcoded params as in old code or reuse given preprocessing options?
  empty = getMeasurementsFD(params.bgParams.emptyMeas, false; bgCorrection = false, numAverages=1, kwargs..., frequencies = frequencies)
  numFrames = acqNumPeriodsPerPatch(params.bgParams.emptyMeas)
  bgFrames = [1+(i-1)*numFrames:i*numFrames for i=1:acqNumPatches(params.bgParams.emptyMeas)]
  for i=1:size(data, 2)  # Is this equivalent to acqNumPatches(bMeas)?
    data[:,i,:] = data[:,i,:] .- mean(empty[:,bgFrames[i],:],dims=(2,3))
  end
  return data
end

function process(algo::MultiPatchReconstructionAlgorithm, params::MultiPatchReconstructionParameter, u::AbstractArray)
  weights = process(algo, params.weightingParams, u, WeightingType(params.weightingParams))

  solver = LeastSquaresParameters(reg = params.reg, solverParams = params.solverParams, weights = weights)

  result = process(algo, solver, algo.ffOp, u)

  return gridresult(result, algo.ffOp.grid, algo.sf)
end

function process(algo::MultiPatchReconstructionAlgorithm, params::Union{W, ProcessResultCache{W}}, u, ::MeasurementBasedWeighting) where W<:AbstractWeightingParameters
  return process(typeof(algo), params, algo.freqs, algo.ffOp, u, algo.arrayType)
end

function process(algo::MultiPatchReconstructionAlgorithm, params::Union{W, ProcessResultCache{W}}, u, ::SystemMatrixBasedWeighting) where W<:AbstractWeightingParameters
  return algo.weights
end