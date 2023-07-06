Base.@kwdef struct SinglePatchReconstructionParameter{L<:AbstractSystemMatrixLoadingParameter, S<:AbstractLinearSolver,
   SP<:AbstractSolverParameters, R<:AbstractRegularization} <: AbstractSinglePatchReconstructionParameters
  # File
  sf::MPIFile
  sfLoad::L
  # Solver
  solver::Type{S}
  solverParams::SP
  reg::Vector{R} = AbstractRegularization[]
  # weightingType::WeightingType = WeightingType.None
end

Base.@kwdef mutable struct SinglePatchReconstructionAlgorithm{P} <: AbstractSinglePatchReconstructionAlgorithm where {P<:AbstractSinglePatchAlgorithmParameters}
  params::P
  # Could also do reconstruction progress meter here
  origParam::Union{AbstractSinglePatchAlgorithmParameters, Nothing} = nothing
  sf::Union{MPIFile, Vector{MPIFile}}
  S::AbstractArray
  grid::RegularGridPositions
  freqs::Vector{Int64}
  output::Channel{Any}
end

function SinglePatchReconstruction(params::SinglePatchParameters{<:AbstractPreProcessingParameters, R, PT}) where {R<:AbstractSinglePatchReconstructionParameters, PT <:AbstractPostProcessingParameters}
  return SinglePatchReconstructionAlgorithm(params)
end
function SinglePatchReconstructionAlgorithm(params::SinglePatchParameters{<:AbstractPreProcessingParameters, R, PT}) where {R<:AbstractSinglePatchReconstructionParameters, PT <:AbstractPostProcessingParameters}
  freqs, S, grid = prepareSystemMatrix(params.reco)
  filter = FrequencyFilteredPreProcessingParameters(freqs, params.pre)
  filteredParams = SinglePatchParameters(filter, params.reco, params.post)
  return SinglePatchReconstructionAlgorithm(filteredParams, params, params.reco.sf, S, grid, freqs, Channel{Any}(Inf))
end
recoAlgorithmTypes(::Type{SinglePatchReconstruction}) = SystemMatrixBasedAlgorithm()
RecoUtils.parameter(algo::SinglePatchReconstructionAlgorithm) = algo.origParam

function prepareSystemMatrix(reco::SinglePatchReconstructionParameter{L,S}) where {L<:AbstractSystemMatrixLoadingParameter, S<:AbstractLinearSolver}
  freqs, sf, grid = process(AbstractMPIReconstructionAlgorithm, reco.sf, reco.sfLoad)
  sf, grid = prepareSF(S, sf, grid) 
  return freqs, sf, grid
end


RecoUtils.take!(algo::SinglePatchReconstructionAlgorithm) = Base.take!(algo.output)

function RecoUtils.put!(algo::SinglePatchReconstructionAlgorithm, data::MPIFile)
  consistenceCheck(algo.sf, data)
  
  result = process(algo, data, algo.params)
  
  # Create Image (maybe image parameter as post params?)
  # TODO make more generic to apply to other pre/reco params as well (pre.numAverage main issue atm)
  pixspacing = (spacing(algo.grid) ./ acqGradient(data)[1] .* acqGradient(algo.sf)[1])*1000u"mm"
  offset = (ffPos(data) .- 0.5 .* calibFov(algo.sf))*1000u"mm" .+ 0.5 .* pixspacing
  dt = acqNumAverages(data)*dfCycle(data)*algo.params.pre.numAverages*1u"s"
  im = makeAxisArray(result, pixspacing, offset, dt)
  result = ImageMeta(im, generateHeaderDict(algo.sf, data))

  Base.put!(algo.output, result)
end

function RecoUtils.similar(algo::SinglePatchReconstructionAlgorithm, data::MPIFile)
  # TODO Check num receive channel

  # This isnt a nice structure atm, because the corrections of the params cant be done individually
  # and in certain parts should only be applied after construction

  # Check if freq still valid
  freqs = algo.freqs
  S = algo.S
  grid = algo.grid
  params = parameter(algo)
  if rxNumSamplingPoints(algo.sf) == rxNumSamplingPoints(data)
    # Ensure that no frequencies are used that are not present in the measurement
    freqParams = fromKwargs(PreProcessedFrequencyFilterParameter; toKwargs([params.pre, params.reco]; flatten = DataType[AbstractSystemMatrixLoadingParameter])...)
    measFreqs = process(AbstractMPIReconstructionAlgorithm, data, freqParams)
    freqs = intersect(algo.freqs, measFreqs)
    if freqs != algo.freqs
      S, grid = getSF(params.reco.sf, freqs, nothing; toKwargs(params.pre)...)
      S, grid = prepareSF(params.reco.solver, S, grid)
    end
  end

  numPeriodGrouping = params.pre.numPeriodGrouping
  if rxNumSamplingPoints(params.reco.sf) > rxNumSamplingPoints(data)
    numPeriodGrouping = rxNumSamplingPoints(params.reco.sf) ÷ rxNumSamplingPoints(data)
  end
  numPeriodAverages = params.pre.numPeriodAverages
  if acqNumPeriodsPerFrame(params.reco.sf) < acqNumPeriodsPerFrame(data)
    numPeriodAverages = acqNumPeriodsPerFrame(data) ÷ (acqNumPeriodsPerFrame(params.reco.sf) * numPeriodGrouping)
  end

  pre = fromKwargs(FrequencyFilteredPreProcessingParameters; toKwargs(params.pre, flatten = DataType[], overwrite = Dict{Symbol, Any}(:numPeriodGrouping => numPeriodGrouping, :numPeriodAverages => numPeriodAverages, :frequencies => freqs))...)
  reco = RecoUtils.similar(algo, data, algo.params.reco)
  post = RecoUtils.similar(algo, data, algo.params.post)

  params = SinglePatchParameters(pre, reco, post)

  result = SinglePatchReconstructionAlgorithm(params, params, algo.params.reco.sf, S, grid, freqs, Channel{Any}(Inf))

  return result
end

function RecoUtils.process(algo::SinglePatchReconstructionAlgorithm, f::MPIFile, params::AbstractPreProcessingParameters)
  result = process(typeof(algo), f, params)
  if eltype(algo.S) != eltype(result)
    @warn "System matrix and measurement have different element data type. Mapping measurment data to system matrix element type."
    result = map(eltype(algo.S),result)
  end
  return result
end


function RecoUtils.process(algo::SinglePatchReconstructionAlgorithm, u::Array, params::SinglePatchReconstructionParameter)
  weights = nothing # getWeights(...)

  B = getLinearOperator(algo, params)

  solver = LeastSquaresParameters(params.solver, B, algo.S, params.reg, params.solverParams)

  result = process(algo, u, solver)

  return gridresult(result, algo.grid, algo.sf)
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter{<:DenseSystemMatixLoadingParameter, S}) where {S}
  return linearOperator(nothing, shape(algo.grid), eltype(algo.S))
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter{<:SparseSystemMatrixLoadingParameter, S}) where {S}
  return linearOperator(params.sfLoad.sparseTrafo, shape(algo.grid), eltype(algo.S))
end