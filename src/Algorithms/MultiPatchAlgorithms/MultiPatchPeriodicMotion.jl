export PeriodicMotionPreProcessing, PeriodicMotionReconstructionParameter
Base.@kwdef struct PeriodicMotionPreProcessing{BG<:AbstractMPIBackgroundCorrectionParameters} <: AbstractMPIPreProcessingParameters{BG}
  # Periodic Motion
  frames::Union{Nothing, UnitRange{Int64}, Vector{Int64}} = nothing
  alpha::Float64 = 3.0
  choosePeak::Int64 = 1
  samplingPrecision::Bool = true
  windowType::Int64 = 1
  higherHarmonic::Int64 = 1
  sf::MultiMPIFile
  tfCorrection::Bool = false
  bgParams::BG = NoBackgroundCorrectionParameters()
  # weightingType::WeightingType = WeightingType.None
end

Base.@kwdef struct PeriodicMotionReconstructionParameter{F<:AbstractFrequencyFilterParameter, S<:AbstractSolverParameters} <: AbstractMultiPatchReconstructionParameters
  sf::MultiMPIFile
  freqFilter::F
  solverParams::S
  λ::Float32
  # weightingType::WeightingType = WeightingType.None
end

function MultiPatchReconstructionAlgorithm(params::MultiPatchParameters{<:PeriodicMotionPreProcessing,<:PeriodicMotionReconstructionParameter,<:AbstractMPIPostProcessingParameters})
  reco = params.reco
  freqs = process(MultiPatchReconstructionAlgorithm, reco.freqFilter, reco.sf)
  return MultiPatchReconstructionAlgorithm(params, nothing, reco.sf, Array, nothing, nothing, nothing, freqs, Channel{Any}(Inf))
end

function AbstractImageReconstruction.put!(algo::MultiPatchReconstructionAlgorithm{MultiPatchParameters{PT, R, T}}, data::MPIFile) where {R, T, PT <: PeriodicMotionPreProcessing}
  result = process(algo, algo.params, data, algo.freqs)

  # Create Image (maybe image parameter as post params?)
  # TODO make more generic to apply to other pre/reco params as well (pre.numAverage main issue atm)
  pixspacing = (voxelSize(algo.sf) ./ sfGradient(data,3) .* sfGradient(algo.sf,3)) * 1000u"mm"
  offset = (fieldOfViewCenter(algo.ffOp.grid) .- 0.5.*fieldOfView(algo.ffOp.grid) .+ 0.5.*spacing(algo.ffOp.grid)) * 1000u"mm"
  dt = acqNumAverages(data) * dfCycle(data) * 1 * 1u"s" # Motion has no averages
  im = makeAxisArray(result, pixspacing, offset, dt)
  result = ImageMeta(im, generateHeaderDict(algo.sf, data))

  Base.put!(algo.output, result)
end

function process(algo::MultiPatchReconstructionAlgorithm, params::PeriodicMotionPreProcessing{NoBackgroundCorrectionParameters},
        f::MPIFile, frequencies::Union{Vector{CartesianIndex{2}}, Nothing} = nothing)
  ffPos_ = ffPos(f)
  motFreq = getMotionFreq(params.sf, f, params.choosePeak) ./ params.higherHarmonic
  tmot = getRepetitionsOfSameState(f, motFreq, params.frames)

  uReco = getMeasurementsMotionCompFD(f, motFreq, tmot, frequencies, params.frames, params.alpha,
    params.samplingPrecision, params.windowType)

  p = numDFPeriodsInMotionCycle(motFreq, params.frames, dfCycle(f))

  mapping = collect(1:acqNumPatches(f)) 
  resortedInd = zeros(Int64, acqNumPatches(f), p)
  
  for i=1:acqNumPatches(f)
    resortedInd[i,:] = unflattenOffsetFieldShift(ffPos_)[i][1:p]
  end

  algo.ffOp = MultiPatchOperator(algo.sf, frequencies,
        #indFFPos=resortedInd[:,1], unused keyword
        FFPos=ffPos_[:,resortedInd[:,1]], mapping=mapping,
        FFPosSF=ffPos_[:,resortedInd[:,1]], bgCorrection = false, tfCorrection = params.tfCorrection)

  return uReco
end

function process(algo::MultiPatchReconstructionAlgorithm,
  params::PeriodicMotionPreProcessing{SimpleExternalBackgroundCorrectionParameters}, f::MPIFile, frequencies::Union{Vector{CartesianIndex{2}}, Nothing} = nothing)
  # Foreground
  fgParams = fromKwargs(PeriodicMotionPreProcessing; toKwargs(params)..., bgParams = NoBackgroundCorrectionParameters())
  result = process(algo, fgParams, f, frequencies)
  # Background
  bgParams = fromKwargs(ExternalPreProcessedBackgroundCorrectionParameters; toKwargs(params)..., bgParams = params.bgParams, spectralLeakageCorrection=true)
  return process(algo, bgParams, result, frequencies)
end

function process(algo::MultiPatchReconstructionAlgorithm, params::PeriodicMotionReconstructionParameter, u::Array)
  solver = LeastSquaresParameters(S = algo.ffOp, reg = [L2Regularization(params.λ)], solverParams = params.solverParams)

  result = process(algo, solver, u)

  return gridresult(result, algo.ffOp.grid, algo.sf)
end