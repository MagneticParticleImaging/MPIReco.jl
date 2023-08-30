export MultiPatchReconstructionAlgorithm, MultiPatchReconstructionParameter
Base.@kwdef struct MultiPatchReconstructionParameter{F<:AbstractFrequencyFilterParameter,O<:AbstractMultiPatchOperatorParameter, S<:AbstractSolverParameters, FF<:AbstractFocusFieldPositions, FFSF<:AbstractFocusFieldPositions} <: AbstractMultiPatchReconstructionParameters
  # File
  sf::MultiMPIFile
  freqFilter::F
  opParams::O
  ffPos::FF = DefaultFocusFieldPositions()
  ffPosSF::FFSF = DefaultFocusFieldPositions()
  solverParams::S
  λ::Float32
  # weightingType::WeightingType = WeightingType.None
end

Base.@kwdef mutable struct MultiPatchReconstructionAlgorithm{P} <: AbstractMultiPatchReconstructionAlgorithm where {P<:AbstractMultiPatchAlgorithmParameters}
  params::P
  # Could also do reconstruction progress meter here
  origParam::Union{AbstractMultiPatchAlgorithmParameters,Nothing} = nothing
  opParams::Union{AbstractMultiPatchOperatorParameter, Nothing} = nothing
  sf::MultiMPIFile
  ffOp::Union{Nothing, MultiPatchOperator}
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
  freqs = process(MultiPatchReconstructionAlgorithm, reco.sf, reco.freqFilter)
  filter = FrequencyFilteredPreProcessingParameters(freqs, params.pre)
  filteredParams = MultiPatchParameters(filter, reco, params.post)

  # Prepare operator construction
  ffPos_ = positions(reco.ffPos)
  ffPosSF = positions(reco.ffPosSF)
  if isnothing(ffPosSF)
    L = length(ffPos(reco.sf[1]))
    ffPosSF = [vec(ffPos(SF))[l] for l=1:L, SF in reco.sf]
  end

  return MultiPatchReconstructionAlgorithm(filteredParams, params, reco.opParams, reco.sf, nothing, ffPos_, ffPosSF, freqs, Channel{Any}(Inf))
end
recoAlgorithmTypes(::Type{MultiPatchReconstruction}) = SystemMatrixBasedAlgorithm()
AbstractImageReconstruction.parameter(algo::MultiPatchReconstructionAlgorithm) = algo.origParam

AbstractImageReconstruction.take!(algo::MultiPatchReconstructionAlgorithm) = Base.take!(algo.output)

function AbstractImageReconstruction.put!(algo::MultiPatchReconstructionAlgorithm, data::MPIFile)
  #consistenceCheck(algo.sf, data)

  algo.ffOp = process(algo, data, algo.opParams)

  result = process(algo, data, algo.params)

  # Create Image (maybe image parameter as post params?)
  # TODO make more generic to apply to other pre/reco params as well (pre.numAverage main issue atm)
  pixspacing = (voxelSize(algo.sf) ./ sfGradient(data,3) .* sfGradient(algo.sf,3)) * 1000u"mm"
  offset = (fieldOfViewCenter(algo.ffOp.grid) .- 0.5.*fieldOfView(algo.ffOp.grid) .+ 0.5.*spacing(algo.ffOp.grid)) * 1000u"mm"
  dt = acqNumAverages(data) * dfCycle(data) * algo.params.pre.pre.numAverages * 1u"s"
  im = makeAxisArray(result, pixspacing, offset, dt)
  result = ImageMeta(im, generateHeaderDict(algo.sf, data))

  Base.put!(algo.output, result)
end

function process(algo::MultiPatchReconstructionAlgorithm, f::MPIFile, params::AbstractMultiPatchOperatorParameter)
  ffPos_ = ffPos(f)
  periodsSortedbyFFPos = unflattenOffsetFieldShift(ffPos_)
  idxFirstPeriod = getindex.(periodsSortedbyFFPos,1)
  ffPos_ = ffPos_[:,idxFirstPeriod]
  gradient = acqGradient(f)[:,:,1,idxFirstPeriod]

  if !isnothing(algo.ffPos)
    ffPos_[:] = algo.ffPos
  end

  ffPosSF = algo.ffPosSF

  return MultiPatchOperator(algo.sf, algo.freqs; toKwargs(params)...,
             gradient = gradient, FFPos = ffPos_, FFPosSF = ffPosSF)
end

function process(t::Type{<:MultiPatchReconstructionAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParameters{NoBackgroundCorrectionParameters, <:CommonPreProcessingParameters})
  kwargs = toKwargs(params, default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))), ignore = [:neglectBGFrames, :bgCorrection])
  result = getMeasurementsFD(f, bgCorrection = false; kwargs...)
  periodsSortedbyFFPos = unflattenOffsetFieldShift(ffPos(f))
  uTotal = Base.similar(result,size(result,1),length(periodsSortedbyFFPos),size(result,3))
  for k=1:length(periodsSortedbyFFPos)
      uTotal[:,k,:] = mean(result[:,periodsSortedbyFFPos[k],:], dims=2)
  end
  return uTotal
end

function process(t::Type{<:MultiPatchReconstructionAlgorithm}, f::MPIFile, params::FrequencyFilteredPreProcessingParameters{SimpleExternalBackgroundCorrectionParameters, <:CommonPreProcessingParameters})
  # Foreground, ignore BGCorrection to reuse preprocessing
  fgParams = CommonPreProcessingParameters(;toKwargs(params)..., bgParams = NoBackgroundCorrectionParameters())
  result = process(t, f, FrequencyFilteredPreProcessingParameters(params.frequencies, fgParams))
  # Background
  kwargs = toKwargs(params, ignore = [:neglectBGFrames, :bgCorrection],
    default = Dict{Symbol, Any}(:frames => params.neglectBGFrames ? (1:acqNumFGFrames(f)) : (1:acqNumFrames(f))))
  bgParams = fromKwargs(FrequencyFilteredBackgroundCorrectionParameters; kwargs..., bgParams = params.bgParams)
  return process(t, result, bgParams)
end

function process(::Type{<:MultiPatchReconstructionAlgorithm}, data::Array, params::FrequencyFilteredBackgroundCorrectionParameters{SimpleExternalBackgroundCorrectionParameters})
  kwargs = toKwargs(params, overwrite = Dict{Symbol, Any}(:frames => params.bgParams.bgFrames), ignore = [:bgParams])
  # TODO migrate with hardcoded params as in old code or reuse given preprocessing options?
  empty = getMeasurementsFD(params.bgParams.emptyMeas, false; bgCorrection = false, numAverages=1, kwargs...)
  numFrames = acqNumPeriodsPerPatch(params.bgParams.emptyMeas)
  bgFrames = [1+(i-1)*numFrames:i*numFrames for i=1:acqNumPatches(params.bgParams.emptyMeas)]
  for i=1:size(data, 2)  # Is this equivalent to acqNumPatches(bMeas)?
    data[:,i,:] = data[:,i,:] .- mean(empty[:,bgFrames[i],:],dims=(2,3))
  end
  return data
end

function process(algo::MultiPatchReconstructionAlgorithm, u::Array, params::MultiPatchReconstructionParameter)
  solver = LeastSquaresParameters(Kaczmarz, nothing, algo.ffOp, [L2Regularization(params.λ)], params.solverParams)

  result = process(algo, u, solver)

  return gridresult(result, algo.ffOp.grid, algo.sf)
end