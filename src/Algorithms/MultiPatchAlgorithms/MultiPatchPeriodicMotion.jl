export PeriodicMotionPreProcessing, PeriodicMotionReconstructionParameter
Base.@kwdef struct PeriodicMotionPreProcessing{BG<:AbstractMPIBackgroundCorrectionParameters, W <: AbstractWeightingParameters} <: AbstractMPIPreProcessingParameters{BG}
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
  weightingParams::Union{W, AbstractUtilityReconstructionParameters{W}} = NoWeightingParameters()
end

Base.@kwdef struct PeriodicMotionReconstructionParameter{F<:AbstractFrequencyFilterParameter, S<:AbstractSolverParameters, R <: AbstractRegularization, arrT <: AbstractArray} <: AbstractMultiPatchReconstructionParameters
  sf::MultiMPIFile
  freqFilter::F
  solverParams::S
  reg::Vector{R} = AbstractRegularization[]
  arrayType::Type{arrT} = Array
end

function MultiPatchReconstructionAlgorithm(params::MultiPatchParameters{<:PeriodicMotionPreProcessing,<:PeriodicMotionReconstructionParameter,<:AbstractMPIPostProcessingParameters})
  reco = params.reco
  freqs = process(MultiPatchReconstructionAlgorithm, reco.freqFilter, reco.sf)
  return MultiPatchReconstructionAlgorithm{typeof(params), reco.arrayType, typeof(reco.arrayType{Float32}(undef, 0))}(params, nothing, reco.sf, nothing, reco.arrayType, nothing, nothing, nothing, freqs, Channel{Any}(Inf))
end

function AbstractImageReconstruction.put!(algo::MultiPatchReconstructionAlgorithm{MultiPatchParameters{PT, R, T}}, data::MPIFile) where {R, T, PT <: PeriodicMotionPreProcessing}
  lock(algo) do 
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
end

function process(algo::MultiPatchReconstructionAlgorithm, params::Union{OP, ProcessResultCache{OP}}, 
  f::MPIFile, frequencies::Union{Vector{CartesianIndex{2}}, Nothing} = nothing) where OP <: PeriodicMotionPreProcessing
  uReco, ffOp, weights = process(typeof(algo), params, f, algo.sf, frequencies)
  algo.ffOp = adapt(algo.arrayType, ffOp)
  algo.weights = adapt(algo.arrayType, weights)
  return adapt(algo.arrayType, uReco)
end

function process(algoT::Type{<:MultiPatchReconstructionAlgorithm}, params::PeriodicMotionPreProcessing{NoBackgroundCorrectionParameters},
        f::MPIFile, sf::MPIFile, frequencies::Union{Vector{CartesianIndex{2}}, Nothing} = nothing)
  @info "Loading Multi Patch motion operator"
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

  # Can't adapt data here because it might be used in background correction
  ffOp = MultiPatchOperator(sf, frequencies,
        #indFFPos=resortedInd[:,1], unused keyword
        FFPos=ffPos_[:,resortedInd[:,1]], mapping=mapping,
        FFPosSF=ffPos_[:,resortedInd[:,1]], bgCorrection = false, tfCorrection = params.tfCorrection)

  weights = process(algoT, params.weightingParams, frequencies, sf, ffOp, nothing, Array)
  return uReco, ffOp, weights
end

function process(algoT::Type{<:MultiPatchReconstructionAlgorithm},
  params::PeriodicMotionPreProcessing{SimpleExternalBackgroundCorrectionParameters}, f::MPIFile, sf::MPIFile, frequencies::Union{Vector{CartesianIndex{2}}, Nothing} = nothing)
  # Foreground
  fgParams = fromKwargs(PeriodicMotionPreProcessing; toKwargs(params)..., bgParams = NoBackgroundCorrectionParameters())
  result, ffOp, weights = process(algoT, fgParams, f, sf, frequencies)
  # Background
  bgParams = fromKwargs(ExternalPreProcessedBackgroundCorrectionParameters; toKwargs(params)..., bgParams = params.bgParams, spectralLeakageCorrection=true)
  return process(algoT, bgParams, result, frequencies), ffOp, weights
end

function process(algo::MultiPatchReconstructionAlgorithm, params::PeriodicMotionReconstructionParameter, u::Array)
  solver = LeastSquaresParameters(S = algo.ffOp, reg = params.reg, solverParams = params.solverParams, weights = algo.weights)

  result = process(algo, solver, u)

  return gridresult(result, algo.ffOp.grid, algo.sf)
end